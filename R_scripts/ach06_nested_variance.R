# ACH-06  Sampling error versus analytical error: nested design. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript ach06_nested_variance.R [data.csv]
# data.csv: site, subsample, reading, value.  Balanced design: same number of subsamples per site and readings per subsample.
set.seed(606)

nested_components <- function(d) {
  d$site <- factor(d$site); d$sub <- factor(paste(d$site, d$subsample)); d$reading <- factor(d$reading)
  a <- nlevels(d$site); b <- nlevels(d$sub) / a; n <- nrow(d) / nlevels(d$sub)
  f <- aov(value ~ site / sub, d); s <- summary(f)[[1]]; ms <- s[["Mean Sq"]]; names(ms) <- trimws(rownames(s))
  ms_site <- ms[1]; ms_sub <- ms[2]; ms_res <- ms[3]
  v_meas <- ms_res; v_sub <- max((ms_sub - ms_res) / n, 0); v_site <- max((ms_site - ms_sub) / (b * n), 0)
  tot <- v_meas + v_sub + v_site
  list(table = s, a = a, b = b, n = n, comp = c(site = v_site, subsample = v_sub, measurement = v_meas), pct = 100 * c(site = v_site, subsample = v_sub, measurement = v_meas) / tot) }
sim_data <- function(sites = 5, subs = 4, reps = 2, sd_site = 30, sd_sub = 12, sd_meas = 3, mu = 250) {
  d <- expand.grid(reading = 1:reps, subsample = 1:subs, site = 1:sites); s <- rnorm(sites, 0, sd_site); u <- rnorm(sites * subs, 0, sd_sub)
  d$value <- round(mu + s[d$site] + u[(d$site - 1) * subs + d$subsample] + rnorm(nrow(d), 0, sd_meas), 1); d[, c("site", "subsample", "reading", "value")] }

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE); d <- if (length(args) >= 1) read.csv(args[1]) else { cat("Practice data: 5 sites, 4 subsamples per site, 2 readings each; true sds 30 (site), 12 (subsample) and 3 (measurement) uS/cm\n\n"); sim_data() }
  r <- nested_components(d); print(r$table)
  cat(sprintf("\nVariance components: site %.1f, subsample %.1f, measurement %.1f\nShare of total variance: site %.0f%%, subsample %.0f%%, measurement %.0f%%\n", r$comp[1], r$comp[2], r$comp[3], r$pct[1], r$pct[2], r$pct[3]))
  cat(sprintf("Standard deviations: site %.1f, subsample %.1f, measurement %.1f\n", sqrt(r$comp[1]), sqrt(r$comp[2]), sqrt(r$comp[3])))
  cat("\nStandard error of the overall mean for different allocations of the same 40 readings (using these components)\n")
  for (cfg in list(c(10, 4, 1), c(5, 4, 2), c(20, 2, 1), c(5, 8, 1))) { a <- cfg[1]; b <- cfg[2]; n <- cfg[3]; se <- sqrt(r$comp[1] / a + r$comp[2] / (a * b) + r$comp[3] / (a * b * n)); cat(sprintf("%2d sites x %d subsamples x %d readings: se = %.2f\n", a, b, n, se)) }
  cat("\nRecovering the true components: 300 simulated experiments, mean estimate of each sd\n")
  est <- t(replicate(300, sqrt(nested_components(sim_data())$comp))); print(round(colMeans(est), 1))
}
