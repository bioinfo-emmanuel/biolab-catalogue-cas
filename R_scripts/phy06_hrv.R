# PHY-06  Heart rate variability from RR intervals. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript phy06_hrv.R [rr.csv]
# rr.csv columns: subject, group, rr_ms  (RR intervals in milliseconds, in time order). With no argument the script
# simulates two groups of 12 people (group B has lower variability) and adds artifacts so it runs offline.
set.seed(66)

clean_rr <- function(rr, lo = 300, hi = 2000, tol = 0.20) {          # range limit, then a 20% rule against the local median
  ok <- rr >= lo & rr <= hi
  med <- as.numeric(stats::filter(rr, rep(1 / 5, 5), sides = 2)); med[is.na(med)] <- median(rr, na.rm = TRUE)
  ok <- ok & abs(rr - med) / med <= tol
  rr[ok]
}
hrv <- function(rr) { d <- diff(rr)
  c(n = length(rr), mean_hr = 60000 / mean(rr), sdnn = sd(rr), rmssd = sqrt(mean(d^2)), pnn50 = 100 * mean(abs(d) > 50)) }

sim_subject <- function(n_beats = 300, mean_rr = 850, rsa = 40, sd_noise = 15, artifacts = 0.02) {
  t <- cumsum(rep(mean_rr, n_beats)) / 1000
  rr <- mean_rr + rsa * sin(2 * pi * 0.25 * t) + rnorm(n_beats, 0, sd_noise)
  k <- sample(n_beats, ceiling(artifacts * n_beats)); rr[k] <- rr[k] * sample(c(0.5, 2), length(k), TRUE)   # missed or extra beat detections
  rr }

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) d <- read.csv(args[1]) else {
    d <- do.call(rbind, lapply(1:24, function(i) { g <- ifelse(i <= 12, "A", "B")
      rr <- sim_subject(mean_rr = rnorm(1, 850, 60), rsa = ifelse(g == "A", 45, 20), sd_noise = ifelse(g == "A", 20, 12)); data.frame(subject = i, group = g, rr_ms = rr) }))
    cat("Practice data: 24 simulated subjects (12 per group), 300 beats each, about 2% artifacts\n") }
  raw <- t(sapply(split(d$rr_ms, d$subject), hrv)); cln <- t(sapply(split(d$rr_ms, d$subject), function(x) hrv(clean_rr(x))))
  grp <- tapply(d$group, d$subject, `[`, 1)
  cat("\nEffect of artifacts on SDNN and RMSSD (medians over subjects)\n")
  print(round(rbind(before_cleaning = c(sdnn = median(raw[, "sdnn"]), rmssd = median(raw[, "rmssd"])), after_cleaning = c(sdnn = median(cln[, "sdnn"]), rmssd = median(cln[, "rmssd"]))), 1))
  cat("\nGroup medians after cleaning\n"); gm <- aggregate(as.data.frame(cln), list(group = grp), median); gm[, -1] <- round(gm[, -1], 1); print(gm, row.names = FALSE)
  for (m in c("sdnn", "rmssd")) { w <- wilcox.test(cln[grp == "A", m], cln[grp == "B", m], exact = FALSE); cat(sprintf("Wilcoxon test for %s, group A vs B: p = %.4f\n", m, w$p.value)) }
  cat(sprintf("Beats removed by the cleaning rule: %.1f%% of all beats\n", 100 * (1 - sum(cln[, "n"]) / sum(raw[, "n"]))))
}
