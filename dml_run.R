# Runs the 20 double machine learning fits of the paper (renal and pulmonary pathways, cross-fitting
# seeds 1001 to 1010), each in its own R process, then pools them as in Results E3: median of the 10
# estimates, with standard error sqrt(median(se^2 + (estimate - median)^2)) (Chernozhukov et al.,
# Econometrics Journal 2018) and a normal 95% interval.
# A fit whose result file already exists is skipped. A fit takes about 8 minutes. The paper's run used
# eight processes in parallel. Each fit sets its own seeds, so the order of the fits does not matter.
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
dir.create("output/dml", showWarnings = FALSE, recursive = TRUE)
for (pathway in c("renal", "pulmonary")) for (seed in 1001:1010) {
  out <- sprintf("output/dml/medoutcon_%s_s%d.csv", pathway, seed)
  if (!file.exists(out)) system2(rscript, c("--vanilla", "dml_fit.R", paste0("pathway=", pathway), paste0("seed=", seed)))
}
files <- list.files("output/dml", pattern = "^medoutcon_.*\\.csv$", full.names = TRUE)
x <- do.call(rbind, lapply(files, utils::read.csv, stringsAsFactors = FALSE))
pooled <- do.call(rbind, lapply(split(x, x[c("pathway", "effect")], drop = TRUE), function(g) {
  med <- stats::median(g$estimate)
  se <- sqrt(stats::median(g$se^2 + (g$estimate - med)^2))
  data.frame(pathway = g$pathway[1], effect = g$effect[1], n_seeds = nrow(g), estimate = med, se = se,
             ci_low = med - stats::qnorm(0.975) * se, ci_high = med + stats::qnorm(0.975) * se,
             min_estimate = min(g$estimate), max_estimate = max(g$estimate))
}))
rownames(pooled) <- NULL
utils::write.csv(pooled, "output/dml_estimates.csv", row.names = FALSE)
print(pooled, digits = 3)
