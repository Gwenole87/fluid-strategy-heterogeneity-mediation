# Double machine learning re-estimate of the two indirect effects at the primary outcome (Results E3),
# with medoutcon (Díaz, Hejazi, Rudolph and van der Laan, Biometrika 2021). One fit is one pathway
# with one cross-fitting seed. dml_run.R runs the 20 fits of the paper and takes the median.
#
# Usage, from the folder that holds the code:
#   Rscript --vanilla dml_fit.R pathway=renal seed=1001
# 10 cross-fitting folds and 5 inner folds are the defaults (folds=, inner=).
# Input: data/mediation/cohort_landmark.csv (renal) or cohort_headroom.csv (pulmonary), written by
# CLOVERS_Mediation.ipynb. Output: output/dml/medoutcon_<pathway>_s<seed>.csv, aggregates only.
#
# Estimands on the risk-difference scale, alive and discharged by day 28, restrictive (1) versus liberal (0):
#   PNDE (pure natural direct effect)    = E[Y(1,M(0))] - E[Y(0,M(0))]
#   TNIE (total natural indirect effect) = E[Y(1,M(1))] - E[Y(1,M(0))]
#   TE   (total effect)                  = E[Y(1,M(1))] - E[Y(0,M(0))]
# Where this departs from the medoutcon tutorials:
#   one call with effect = "pm" returns theta(1,1), theta(1,0) and theta(0,0) with their efficient
#   influence functions, and the three contrasts are formed from them,
#   Z = NULL (no intermediate confounder) gives the natural effects,
#   sl3's own xgboost learner calls functions removed in xgboost 3, hence Lrnr_xgb_es below,
#   the fit is repeated over 10 cross-fitting seeds and the median is reported (dml_run.R).
# Patient-level data are imputed in memory and never written to disk.

Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")

DATA_DIR <- "data/mediation"
OUT_DIR  <- "output/dml"

EXPOSURE <- "restrictive"
OUTCOME  <- "alive_discharged_28d"
BASEC <- c("age", "female", "sofa_baseline", "lactate_baseline", "creatinine_baseline", "sf_ratio_baseline",
           "map_baseline", "ckd_any", "sepsis_source_pulmonary", "on_vasopressor_baseline", "prerand_fluid_L")
PATHWAYS <- list(
  renal     = list(file = "cohort_landmark.csv", n = 1504L, mediator = "delta_creat_0_24h",   covars = BASEC),
  pulmonary = list(file = "cohort_headroom.csv", n = 1408L, mediator = "resp_escalation_24h", covars = c(BASEC, "on_intermediate_baseline"))
)
IMPUTATION_SEED <- 20261027L

# Learner library: main-terms GLM and xgboost at depths 1, 2, 4 and 6, each xgboost with early stopping
# on the number of trees (20% internal hold-out of its own training data, patience 20 rounds, at most
# 1,000 trees, prediction with the best iteration). For each nuisance model, the learner with the best
# 5-fold cross-validated loss inside the cross-fitting training set is selected (discrete super learner).
LEARNER_SETTINGS <- list(depths = c(1L, 2L, 4L, 6L), eta = 0.05, nrounds_max = 1000L, early_stopping_rounds = 20L,
                         valid_frac = 0.2, subsample = 0.8, colsample_bytree = 0.8, min_child_weight = 1, inner_folds = 5L)

# xgboost 3: best_iteration is 0-based, so the number of rounds kept for prediction is best_iteration + 1.
XGB_BEST_ROUNDS <- function(fit) as.integer(xgboost::xgb.attr(fit, "best_iteration")) + 1L

load_cohort <- function(pathway) {
  p <- PATHWAYS[[pathway]]
  df <- utils::read.csv(file.path(DATA_DIR, p$file), stringsAsFactors = FALSE)
  stopifnot(nrow(df) == p$n, all(c(EXPOSURE, OUTCOME, p$mediator, p$covars) %in% names(df)))
  df <- df[, c(EXPOSURE, p$covars, p$mediator, OUTCOME)]
  stopifnot(!anyNA(df[[EXPOSURE]]), !anyNA(df[[OUTCOME]]), all(df[[EXPOSURE]] %in% 0:1), all(df[[OUTCOME]] %in% 0:1))
  df
}

# Single imputation with missingness indicators. Baseline covariates are imputed once (predictive mean
# matching) from the other baseline covariates only, with a 0/1 indicator for every covariate that has a
# missing value. The mediator is imputed once from strategy, covariates and outcome (predictive mean
# matching if continuous, logistic regression if binary).
impute_cohort <- function(df, pathway, m = 1L, seed = IMPUTATION_SEED, maxit = 10L) {
  p <- PATHWAYS[[pathway]]
  cov_miss <- p$covars[colSums(is.na(df[, p$covars, drop = FALSE])) > 0]
  ind <- if (length(cov_miss)) stats::setNames(as.data.frame(lapply(cov_miss, function(v) as.integer(is.na(df[[v]])))),
                                               paste0("miss_", cov_miss)) else NULL
  med <- p$mediator
  med_missing <- anyNA(df[[med]])
  med_binary <- all(stats::na.omit(df[[med]]) %in% 0:1)
  if (!length(cov_miss) && !med_missing) return(rep(list(if (is.null(ind)) df else cbind(df, ind)), m))
  vars <- names(df)
  pm <- matrix(0L, length(vars), length(vars), dimnames = list(vars, vars))
  meth <- stats::setNames(rep("", length(vars)), vars)
  for (v in cov_miss) { pm[v, setdiff(p$covars, v)] <- 1L; meth[v] <- "pmm" }
  dfi <- df
  if (med_missing) {
    pm[med, c(EXPOSURE, p$covars, OUTCOME)] <- 1L
    meth[med] <- if (med_binary) "logreg" else "pmm"
    if (med_binary) dfi[[med]] <- factor(dfi[[med]], levels = c(0, 1))
  }
  imp <- mice::mice(dfi, m = m, method = meth, predictorMatrix = pm, maxit = maxit, seed = seed, printFlag = FALSE)
  lapply(seq_len(m), function(i) {
    d <- mice::complete(imp, i)
    if (med_missing && med_binary) d[[med]] <- as.integer(as.character(d[[med]]))
    if (!is.null(ind)) d <- cbind(d, ind)
    stopifnot(!anyNA(d), nrow(d) == nrow(df))
    d
  })
}

analysis_columns <- function(d, pathway) {
  p <- PATHWAYS[[pathway]]
  list(W = c(p$covars, grep("^miss_", names(d), value = TRUE)), A = EXPOSURE, M = p$mediator, Y = OUTCOME)
}

# sl3 learners for medoutcon
make_Lrnr_xgb_es <- function() {
  R6::R6Class(
    classname = "Lrnr_xgb_es", inherit = sl3::Lrnr_base, portable = TRUE, class = TRUE,
    public = list(
      initialize = function(max_depth = 2L, eta = 0.05, nrounds = 1000L, early_stopping_rounds = 20L, valid_frac = 0.2,
                            subsample = 0.8, colsample_bytree = 0.8, min_child_weight = 1, nthread = 1L, ...) {
        params <- list(max_depth = max_depth, eta = eta, nrounds = nrounds, early_stopping_rounds = early_stopping_rounds,
                       valid_frac = valid_frac, subsample = subsample, colsample_bytree = colsample_bytree,
                       min_child_weight = min_child_weight, nthread = nthread, ...)
        super$initialize(params = params, ...)
      }
    ),
    private = list(
      .properties = c("continuous", "binomial", "weights"),
      .train = function(task) {
        p <- self$params
        outcome_type <- self$get_outcome_type(task)
        binary <- identical(outcome_type$type, "binomial")
        y <- as.numeric(outcome_type$format(task$Y))
        X <- as.matrix(task$X); storage.mode(X) <- "double"
        w <- if (task$has_node("weights")) task$weights else rep(1, length(y))
        n <- length(y)
        va <- if (binary) {
          unlist(lapply(split(seq_len(n), y), function(ix) ix[sample.int(length(ix), max(1L, round(p$valid_frac * length(ix))))]))
        } else sample.int(n, max(1L, round(p$valid_frac * n)))
        tr <- setdiff(seq_len(n), va)
        prm <- list(objective = if (binary) "binary:logistic" else "reg:squarederror",
                    eval_metric = if (binary) "logloss" else "rmse",
                    eta = p$eta, max_depth = p$max_depth, subsample = p$subsample,
                    colsample_bytree = p$colsample_bytree, min_child_weight = p$min_child_weight, nthread = p$nthread)
        dtr <- xgboost::xgb.DMatrix(X[tr, , drop = FALSE], label = y[tr], weight = w[tr], nthread = p$nthread)
        dva <- xgboost::xgb.DMatrix(X[va, , drop = FALSE], label = y[va], weight = w[va], nthread = p$nthread)
        fit <- xgboost::xgb.train(params = prm, data = dtr, nrounds = p$nrounds, evals = list(valid = dva),
                                  early_stopping_rounds = p$early_stopping_rounds, verbose = 0)
        best <- XGB_BEST_ROUNDS(fit)
        list(model = fit, best = best, feature_names = colnames(X))
      },
      .predict = function(task) {
        fo <- private$.fit_object
        X <- as.matrix(task$X)[, fo$feature_names, drop = FALSE]; storage.mode(X) <- "double"
        as.numeric(stats::predict(fo$model, xgboost::xgb.DMatrix(X, nthread = 1L), iterationrange = c(1L, fo$best)))
      },
      .required_packages = c("xgboost")
    )
  )
}

make_sl3_library <- function(s = LEARNER_SETTINGS) {
  Lrnr_xgb_es <- make_Lrnr_xgb_es()
  base <- c(list(sl3::Lrnr_glm$new()),
            lapply(s$depths, function(d) Lrnr_xgb_es$new(max_depth = d, eta = s$eta, nrounds = s$nrounds_max,
                   early_stopping_rounds = s$early_stopping_rounds, valid_frac = s$valid_frac, subsample = s$subsample,
                   colsample_bytree = s$colsample_bytree, min_child_weight = s$min_child_weight, nthread = 1L)))
  stack <- do.call(sl3::Stack$new, base)
  list(
    bin  = sl3::Lrnr_sl$new(learners = stack, cv_control = list(V = s$inner_folds),
                            metalearner = sl3::Lrnr_cv_selector$new(eval_function = sl3::loss_loglik_binomial)),
    cont = sl3::Lrnr_sl$new(learners = stack, cv_control = list(V = s$inner_folds),
                            metalearner = sl3::Lrnr_cv_selector$new(eval_function = sl3::loss_squared_error))
  )
}

effects_from_thetas <- function(th, eif, n) {
  # th: named c(t11, t10, t00); eif: list of centred EIF vectors with the same names
  est <- c(PNDE = th[["t10"]] - th[["t00"]], TNIE = th[["t11"]] - th[["t10"]], TE = th[["t11"]] - th[["t00"]])
  se  <- c(PNDE = stats::sd(eif$t10 - eif$t00), TNIE = stats::sd(eif$t11 - eif$t10), TE = stats::sd(eif$t11 - eif$t00)) / sqrt(n)
  data.frame(effect = names(est), estimate = unname(est), se = unname(se))
}

fit_medoutcon_once <- function(d, cols, K, seed, s = LEARNER_SETTINGS, estimator = "onestep") {
  lib <- make_sl3_library(s)
  set.seed(seed)
  fit <- medoutcon::medoutcon(
    W = as.data.frame(d[, cols$W]), A = d[[cols$A]], Z = NULL, M = as.data.frame(d[, cols$M, drop = FALSE]), Y = d[[cols$Y]],
    effect = "pm",                                  # returns theta(1,1), theta(0,0), theta(1,0) with their EIFs
    g_learners = lib$bin, h_learners = lib$bin, b_learners = lib$bin, v_learners = lib$cont,
    u_learners = sl3::Lrnr_mean$new(),              # u term is multiplied by (Z - q) = 0 when Z = NULL
    d_learners = sl3::Lrnr_mean$new(),              # two-phase sampling term, not used (all mediators observed)
    estimator = estimator, estimator_args = list(cv_folds = K, cv_strat = FALSE, strat_pmin = 0.1))
  cr <- fit$contrast_results
  th  <- c(t11 = cr$contrast_1_1_mean, t10 = cr$contrast_1_0_mean, t00 = cr$contrast_0_0_mean)
  eif <- list(t11 = cr$contrast_1_1_eif, t10 = cr$contrast_1_0_eif, t00 = cr$contrast_0_0_eif)
  out <- effects_from_thetas(th, eif, nrow(d))
  rbind(out, data.frame(effect = c("E[Y(1,M(1))]", "E[Y(1,M(0))]", "E[Y(0,M(0))]"), estimate = unname(th),
                        se = sapply(eif, stats::sd) / sqrt(nrow(d))))
}

args <- commandArgs(trailingOnly = TRUE)
kv <- function(key, default = NULL) {
  hit <- grep(paste0("^", key, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^", key, "="), "", hit[1]) else default
}
pathway <- kv("pathway", "renal"); stopifnot(pathway %in% names(PATHWAYS))
seed <- as.integer(kv("seed", "1001")); K <- as.integer(kv("folds", "10")); inner <- as.integer(kv("inner", "5"))
s <- LEARNER_SETTINGS; s$inner_folds <- inner
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages({
  library(data.table); data.table::setDTthreads(1L)
  library(sl3); library(medoutcon)
})
df <- load_cohort(pathway)
d <- impute_cohort(df, pathway, m = 1L)[[1]]
cols <- analysis_columns(d, pathway)
res <- fit_medoutcon_once(d, cols, K, seed, s, estimator = "onestep")
res <- cbind(data.frame(package = "medoutcon", pathway = pathway, folds = K, inner_folds = inner, seed = seed,
                        n = nrow(d), stringsAsFactors = FALSE), res)
utils::write.csv(res, file.path(OUT_DIR, sprintf("medoutcon_%s_s%d.csv", pathway, seed)), row.names = FALSE)
cat(pathway, "seed", seed, "done\n")
