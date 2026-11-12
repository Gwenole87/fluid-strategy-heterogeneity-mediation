# Odds ratios, risk ratios and risk differences from the counterfactual probabilities saved by
# cmaverse_models.R (output/cmaverse/model_XX.rds), with the pooling of CMAverse 0.1.0: the point
# estimate is the mean over the 10 imputations of the original data (on the log scale for
# ratios), the 95% interval is given by the 2.5th and 97.5th percentiles of the 1,000 bootstrap
# replicates (each the mean over its 10 imputations), and p = 2 x min(share above, share below
# the null). E-values on the risk-ratio scale (VanderWeele and Ding), for the estimate and for
# the confidence limit closest to 1.
suppressPackageStartupMessages(library(EValue))
RD <- "output/cmaverse"
files <- sort(list.files(RD, pattern = "^model_[0-9]{2}\\.rds$", full.names = TRUE))
lor <- function(p1, p0) log(p1 / (1 - p1)) - log(p0 / (1 - p0))
lrr <- function(p1, p0) log(p1) - log(p0)
rdf <- function(p1, p0) p1 - p0
contr <- list(cde = c("EY1m", "EY0m"), pnde = c("EY10", "EY00"), tnde = c("EY11", "EY01"),
              pnie = c("EY01", "EY00"), tnie = c("EY11", "EY10"), te = c("EY11", "EY00"))
out <- list(); means <- list()
for (f in files) {
  m <- readRDS(f); cap <- m$captured_means
  stopifnot(nrow(cap) == 10L + 10L + 10000L)
  pt <- cap[1:10, , drop = FALSE]; bt <- cap[21:nrow(cap), , drop = FALSE]
  rep_id <- rep(seq_len(1000), each = 10)
  means[[length(means) + 1]] <- data.frame(model = m$k, label = m$label, t(colMeans(pt)))
  for (en in names(contr)) {
    a <- contr[[en]][1]; b0 <- contr[[en]][2]
    for (sc in c("OR", "RR", "RD")) {
      fn <- switch(sc, OR = lor, RR = lrr, RD = rdf)
      est_imp <- fn(pt[, a], pt[, b0]); est <- mean(est_imp)
      reps <- tapply(fn(bt[, a], bt[, b0]), rep_id, mean)
      ci <- quantile(reps, c(0.025, 0.975), na.rm = TRUE)
      pv <- if (est == 0) 1 else 2 * min(sum(reps > 0), sum(reps < 0)) / sum(!is.na(reps))
      tr <- if (sc == "RD") identity else exp
      out[[length(out) + 1]] <- data.frame(model = m$k, label = m$label, effect = en, scale = sc,
                                           estimate = tr(est), lo = tr(ci[[1]]), hi = tr(ci[[2]]), p = pv)
    }
  }
}
res <- do.call(rbind, out); mm <- do.call(rbind, means)
# E-values on the risk-ratio scale, and on the odds-ratio scale with the conversion for a common outcome (rare = FALSE)
ev <- function(est, lo, hi, sc) {
  if (sc == "RR") e <- evalues.RR(est = est, lo = lo, hi = hi)
  else if (sc == "OR") e <- evalues.OR(est = est, lo = lo, hi = hi, rare = FALSE)
  else return(c(NA, NA))
  c(e["E-values", "point"], suppressWarnings(min(e["E-values", c("lower", "upper")], na.rm = TRUE)))
}
evs <- t(mapply(ev, res$estimate, res$lo, res$hi, res$scale))
res$evalue_point <- evs[, 1]; res$evalue_ci <- evs[, 2]
write.csv(res, file.path(RD, "effect_scales.csv"), row.names = FALSE)
write.csv(mm, file.path(RD, "counterfactual_means.csv"), row.names = FALSE)
cat("models processed:", length(files), "\n")
