#!/usr/bin/env Rscript
# CLOVERS 24-hour landmark mediation: CMAverse g-formula, natural direct and indirect effects.
#
# This is the script that CLOVERS_Mediation.ipynb writes and runs (section D2), with one addition.
# At every g-formula evaluation, the six counterfactual probabilities computed by CMAverse
# (EY0m, EY1m, EY00, EY01, EY10, EY11) are saved: for the original data and for each of the
# 1,000 bootstrap resamples, in each of the 10 imputations. effect_scales.R turns them into risk
# ratios, risk differences and odds ratios. The addition draws no random number, so the odds
# ratios are the same as in the notebook run.
#
# Input: data/mediation/cohort_headroom.csv and cohort_landmark.csv, written by the notebook.
# Output: output/cmaverse/model_01.rds to model_14.rds, and cmaverse_results.json.
# The 14 models take about 18 hours. A checkpoint is saved after each model, and a stopped run
# restarts where it stopped, with the same random-number state.
#
# Joint two-mediator model, single-mediator models and alternative analyses. The pulmonary models
# use the escalation cohort (called headroom here: baseline respiratory level 0 or 1) and adjust
# for on_intermediate_baseline (BASEC_PULM). postc is empty: no post-randomisation confounder is
# adjusted for.

suppressPackageStartupMessages({ library(CMAverse); library(jsonlite) })
# Capture of the six counterfactual probabilities at every est.gformula call (no random number
# is drawn), and a checkpoint after every model so that the sequence can be resumed.
.cap_env <- new.env(); .cap_env$log <- list()
local({
  ns <- asNamespace("CMAverse"); f <- get("est.gformula", envir = ns); b <- as.list(body(f))
  idx <- which(sapply(b, function(e) is.call(e) && identical(e[[1]], as.name("rm")) &&
                 any(sapply(as.list(e)[-1], function(z) identical(z, as.name("weightsEY"))))))
  stopifnot(length(idx) == 1)
  cap <- quote({ .cap_env$log[[length(.cap_env$log) + 1L]] <- c(EY0m = EY0m, EY1m = EY1m, EY00 = EY00, EY01 = EY01, EY10 = EY10, EY11 = EY11)
    if (length(.cap_env$log) == 10L) saveRDS(do.call(rbind, .cap_env$log), file.path(OUT_DIR, sprintf("model_%02d_point.rds", .MODEL_K))) })
  body(f) <- as.call(append(b, list(cap), after = idx - 1))
  assignInNamespace("est.gformula", f, ns = "CMAverse")
})
OUT_DIR <- "output/cmaverse"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
CKPT <- file.path(OUT_DIR, "checkpoint.rds")
.DONE <- 0L; .SAVED <- list(); .SEED_RESUME <- NULL; .MODEL_K <- 0L
if (file.exists(CKPT)) { ck <- readRDS(CKPT); .DONE <- ck$done; .SAVED <- ck$saved; .SEED_RESUME <- ck$seed
  cat(sprintf("checkpoint found: %d models done\n", .DONE)) }
save_model <- function(k, label, res, eff, ev, keep, t_start) {
  cap <- do.call(rbind, .cap_env$log)
  el <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  saveRDS(list(k = k, label = label, basec = keep, n = if (is.null(res)) NA else nrow(res$data),
               effect.pe = if (is.null(res)) NULL else res$effect.pe, effect.ci.low = if (is.null(res)) NULL else res$effect.ci.low,
               effect.ci.high = if (is.null(res)) NULL else res$effect.ci.high, effect.pval = if (is.null(res)) NULL else res$effect.pval,
               evalues = ev, captured_means = cap, elapsed_s = el, finished = format(Sys.time())),
          file.path(OUT_DIR, sprintf("model_%02d.rds", k)))
  ck <- list(done = k, saved = RESULTS, seed = get(".Random.seed", envir = globalenv()))
  saveRDS(ck, CKPT); saveRDS(ck, file.path(OUT_DIR, sprintf("checkpoint_after_model_%02d.rds", k)))
  cat(sprintf("  model %d saved, %d captured calls, %.0f s\n", k, if (is.null(cap)) 0L else nrow(cap), el))
}
cat(sprintf("CMAverse %s\n", packageVersion("CMAverse")))

set.seed(42)
NBOOT <- 1000L
NREP  <- 100L    # passed to cmest, not used by the g-formula estimator
NMICE <- 10L    # mice imputations for missing data

df_headroom <- read.csv("data/mediation/cohort_headroom.csv", stringsAsFactors = FALSE)
df_landmark <- read.csv("data/mediation/cohort_landmark.csv", stringsAsFactors = FALSE)
df_headroom$restrictive <- as.integer(df_headroom$restrictive)
df_landmark$restrictive <- as.integer(df_landmark$restrictive)
cat(sprintf("headroom n=%d | landmark n=%d\n", nrow(df_headroom), nrow(df_landmark)))

BASEC      <- c("age", "female", "sofa_baseline", "lactate_baseline", "creatinine_baseline", "sf_ratio_baseline", "map_baseline", "ckd_any", "sepsis_source_pulmonary", "on_vasopressor_baseline", "prerand_fluid_L")
BASEC_PULM <- c("age", "female", "sofa_baseline", "lactate_baseline", "creatinine_baseline", "sf_ratio_baseline", "map_baseline", "ckd_any", "sepsis_source_pulmonary", "on_vasopressor_baseline", "prerand_fluid_L", "on_intermediate_baseline")
BASEC_ROB  <- c("age", "female", "sofa_baseline", "lactate_baseline", "creatinine_baseline", "sf_ratio_baseline", "map_baseline", "ckd_any", "sepsis_source_pulmonary", "on_vasopressor_baseline", "prerand_fluid_L", "on_adv_support_baseline")

# Helpers
extract_effects <- function(res) {
  pe <- res$effect.pe; cl <- res$effect.ci.low; ch <- res$effect.ci.high; pv <- res$effect.pval
  if (is.null(pe) || length(pe) == 0) return(list())
  eff <- list()
  for (en in names(pe)) {
    eff[[en]] <- list(estimate = as.numeric(pe[[en]]),
                      ci_lower = as.numeric(cl[[en]]),
                      ci_upper = as.numeric(ch[[en]]),
                      p_value  = if (!is.null(pv) && en %in% names(pv)) as.numeric(pv[[en]]) else NA_real_)
  }
  eff
}

extract_evalues <- function(res) {
  tryCatch({
    sens <- tryCatch(cmsens(object = res, sens = "uc"),
                     error = function(e1) tryCatch(cmsens(object = res, sens.type = "uc"),
                                                   error = function(e2) NULL))
    if (is.null(sens) || is.null(sens$evalues)) return(list())
    ev_tbl <- sens$evalues
    evcols <- grep("^Evalue", colnames(ev_tbl), value = TRUE)
    ev_est_col <- evcols[grep("est",   evcols)[1]]
    ev_ci_col  <- evcols[grep("lower", evcols)[1]]
    if (is.na(ev_est_col) && ncol(ev_tbl) >= 4L) ev_est_col <- colnames(ev_tbl)[4]
    if (is.na(ev_ci_col)  && ncol(ev_tbl) >= 5L) ev_ci_col  <- colnames(ev_tbl)[5]
    evals <- list()
    for (rn in rownames(ev_tbl)) {
      est <- if (!is.na(ev_est_col)) as.numeric(ev_tbl[rn, ev_est_col]) else NA_real_
      ci  <- if (!is.na(ev_ci_col))  as.numeric(ev_tbl[rn, ev_ci_col])  else NA_real_
      if (!is.na(est)) evals[[rn]] <- list(evalue_est = est, evalue_ci = ci)
    }
    evals
  }, error = function(e) { cat(sprintf("    E-value error: %s\n", e$message)); list() })
}

RESULTS <- list()

run_cmest <- function(df, mediators, mregs, outcome, yreg, basec_use, label) {
  cat("\n"); cat(strrep("=", 70)); cat("\n"); cat(label); cat("\n"); cat(strrep("=", 70)); cat("\n")
  .MODEL_K <<- .MODEL_K + 1L
  if (.MODEL_K <= .DONE) { RESULTS[[length(RESULTS) + 1]] <<- .SAVED[[.MODEL_K]]; cat("  skipped, done in checkpoint\n"); return(invisible(NULL)) }
  if (.MODEL_K == .DONE + 1L && .DONE > 0L) { assign(".Random.seed", .SEED_RESUME, envir = globalenv()); cat("  RNG state restored from checkpoint\n") }
  .cap_env$log <- list(); .t_start <- Sys.time()
  if (!all(c(mediators, outcome) %in% colnames(df))) { cat("  SKIP: missing columns\n"); return(invisible(NULL)) }
  basec_use <- basec_use[basec_use %in% colnames(df)]
  keep <- basec_use[sapply(basec_use, function(cc) length(unique(df[[cc]][!is.na(df[[cc]])])) > 1)]
  dropped <- setdiff(basec_use, keep)
  if (length(dropped) > 0) cat(sprintf("  pruned constant basec: %s\n", paste(dropped, collapse = ", ")))
  mcols <- c("restrictive", keep, mediators, outcome); mcols <- mcols[mcols %in% colnames(df)]
  max_miss <- max(colMeans(is.na(df[, mcols, drop = FALSE])))
  use_mi <- max_miss > 0.05
  cat(sprintf("  n=%d  max missing=%.1f%%  MI=%s (mice m=%d)  nRep=%d  nboot=%d\n",
              nrow(df), max_miss * 100, ifelse(use_mi, "YES", "NO"),
              ifelse(use_mi, NMICE, 0L), NREP, NBOOT))
  mval <- lapply(seq_along(mediators), function(i) {
    if (mregs[[i]] == "logistic") 0L else as.numeric(median(df[[mediators[i]]], na.rm = TRUE))
  })
  args <- list(data = df, model = "gformula", outcome = outcome, exposure = "restrictive",
               mediator = mediators, basec = keep, EMint = TRUE, mreg = mregs, yreg = yreg,
               astar = 0L, a = 1L, mval = mval,
               estimation = "imputation", inference = "bootstrap", nboot = NBOOT, nRep = NREP)
  if (use_mi) { args$multimp <- TRUE; args$args_mice <- list(m = NMICE) }
  res <- tryCatch(do.call(cmest, args),
                  error = function(e) { cat(sprintf("  FAILED: %s\n", e$message)); NULL })
  if (is.null(res)) { RESULTS[[length(RESULTS) + 1]] <<- list(label = label, error = TRUE); save_model(.MODEL_K, label, NULL, NULL, NULL, keep, .t_start); return(invisible(NULL)) }
  eff <- extract_effects(res); ev <- extract_evalues(res)
  for (en in grep("^R", names(eff), value = TRUE)) { e <- eff[[en]]   # ratio effects only
    cat(sprintf("  %-12s %.4f (%.4f, %.4f) p=%.4f\n",
                en, e$estimate, e$ci_lower, e$ci_upper, ifelse(is.na(e$p_value), NA, e$p_value))) }
  RESULTS[[length(RESULTS) + 1]] <<- list(label = label, mediators = as.list(mediators),
      outcome = outcome, yreg = yreg, n = nrow(res$data), basec = as.list(keep),
      effects = eff, evalues = ev)
  save_model(.MODEL_K, label, res, eff, ev, keep, .t_start)
  invisible(res)
}

# Primary outcome
# Joint two-mediator natural-effects model (headroom): PNDE + joint indirect
run_cmest(df_headroom, c("delta_creat_0_24h", "resp_escalation_24h"), list("linear", "logistic"),
          "alive_discharged_28d", "logistic", BASEC_PULM,
          "JOINT two-mediator | 28d alive-discharged [joint: PNDE + joint indirect]")

# Single renal (full landmark): indirect via renal
run_cmest(df_landmark, c("delta_creat_0_24h"), list("linear"),
          "alive_discharged_28d", "logistic", BASEC,
          "SINGLE renal (delta-creatinine) | 28d alive-discharged [indirect-via-renal]")

# Single pulmonary (headroom, gradient escalation): indirect via pulmonary
run_cmest(df_headroom, c("resp_escalation_24h"), list("logistic"),
          "alive_discharged_28d", "logistic", BASEC_PULM,
          "SINGLE pulmonary (gradient escalation) | 28d alive-discharged [indirect-via-pulmonary]")

# Single renal on the headroom cohort (compared with the joint model in the notebook, Part E)
run_cmest(df_headroom, c("delta_creat_0_24h"), list("linear"),
          "alive_discharged_28d", "logistic", BASEC_PULM,
          "SINGLE renal (headroom) | 28d alive-discharged [overlap-cohort]")

# Secondary outcome (90-day mortality)
run_cmest(df_headroom, c("delta_creat_0_24h", "resp_escalation_24h"), list("linear", "logistic"),
          "mortality_90d", "logistic", BASEC_PULM,
          "JOINT two-mediator | 90d mortality [secondary]")
run_cmest(df_landmark, c("delta_creat_0_24h"), list("linear"),
          "mortality_90d", "logistic", BASEC,
          "SINGLE renal | 90d mortality [secondary]")
run_cmest(df_headroom, c("resp_escalation_24h"), list("logistic"),
          "mortality_90d", "logistic", BASEC_PULM,
          "SINGLE pulmonary | 90d mortality [secondary]")

# Alternative analyses (Results E1 and E2)
# Pulmonary: new advanced support among patients at level 0 at baseline
run_cmest(subset(df_headroom, base_resp_lvl == 0), c("new_adv_support_24h"), list("logistic"),
          "alive_discharged_28d", "logistic", BASEC,
          "ALT pulmonary de-novo composite (level-0 only) | 28d alive-discharged")
# Pulmonary: IMV-only nested rung, headroom cohort
run_cmest(df_headroom, c("new_imv_24h"), list("logistic"),
          "alive_discharged_28d", "logistic", BASEC_PULM,
          "ALT pulmonary IMV-only (nested rung) | 28d alive-discharged")
# Pulmonary: change in SpO2/FiO2 as mediator
run_cmest(df_headroom, c("delta_sf_ratio_0_24h"), list("linear"),
          "alive_discharged_28d", "logistic", BASEC_PULM,
          "ALT pulmonary delta-SpO2/FiO2 (ceiling) | 28d alive-discharged")
# Renal: on_adv_support_baseline added to basec (full landmark)
run_cmest(df_landmark, c("delta_creat_0_24h"), list("linear"),
          "alive_discharged_28d", "logistic", BASEC_ROB,
          "ALT renal + on_adv_support_baseline in basec | 28d alive-discharged")
# Renal pathway without chronic dialysis, without early renal replacement therapy, or without both
run_cmest(subset(df_landmark, chronic_dialysis == 0), c("delta_creat_0_24h"), list("linear"),
          "alive_discharged_28d", "logistic", BASEC,
          "ALT renal excl. baseline dialysis | 28d alive-discharged")
run_cmest(subset(df_landmark, rrt_0_24h == 0), c("delta_creat_0_24h"), list("linear"),
          "alive_discharged_28d", "logistic", BASEC,
          "ALT renal excl. incident RRT 0-24h | 28d alive-discharged")
run_cmest(subset(df_landmark, chronic_dialysis == 0 & rrt_0_24h == 0), c("delta_creat_0_24h"), list("linear"),
          "alive_discharged_28d", "logistic", BASEC,
          "ALT renal excl. baseline dialysis AND incident RRT | 28d alive-discharged")

# Save
writeLines(toJSON(RESULTS, auto_unbox = TRUE, pretty = TRUE), file.path(OUT_DIR, "cmaverse_results.json"))
cat(sprintf("\nSaved %d model results: %s\n", length(RESULTS), file.path(OUT_DIR, "cmaverse_results.json")))
