# ENV-04  Earthquake magnitude-frequency statistics and the b-value. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Input: a CSV with a column named mag (USGS catalog exports use this name).
# With no file given, the script uses R's built-in `quakes` data (1000 events near Fiji, 1964 onward).
args <- commandArgs(trailingOnly = TRUE)
mags <- if (length(args) >= 1) read.csv(args[1])$mag else datasets::quakes$mag
mags <- mags[!is.na(mags)]

# Magnitude of completeness by the maximum-curvature method: the modal bin of the frequency-magnitude curve.
mc_maxcurv <- function(m, dm = 0.1) { b <- round(m / dm) * dm; t <- table(b); as.numeric(names(t)[which.max(t)]) }

# Aki (1965) maximum-likelihood b-value with Utsu's correction for magnitudes binned at dm.
b_mle <- function(m, mc, dm = 0.1) { m <- m[m >= mc]; n <- length(m)
  b <- log10(exp(1)) / (mean(m) - (mc - dm / 2)); c(b = b, se = b / sqrt(n), n = n) }

# Least-squares fit to the cumulative curve, shown for comparison (known to be biased).
b_lsq <- function(m, mc) { m <- m[m >= mc]; x <- sort(unique(m)); N <- sapply(x, function(v) sum(m >= v))
  unname(-coef(lm(log10(N) ~ x))[2]) }

# Bootstrap interval for b.
b_boot <- function(m, mc, B = 1000, dm = 0.1) { m <- m[m >= mc]
  quantile(replicate(B, b_mle(sample(m, replace = TRUE), mc, dm)["b"]), c(0.025, 0.975)) }

if (sys.nframe() == 0) {
  set.seed(3)
  mc <- mc_maxcurv(mags)
  cat(sprintf("Events: %d   magnitude range %.1f to %.1f\nMc (maximum curvature): %.1f\n", length(mags), min(mags), max(mags), mc))
  r <- b_mle(mags, mc); cat(sprintf("b (maximum likelihood) = %.2f +/- %.2f from %d events above Mc\n", r["b"], r["se"], r["n"]))
  cat(sprintf("b (least squares, cumulative) = %.2f\n", b_lsq(mags, mc)))
  cat("95% bootstrap interval for b:", round(b_boot(mags, mc), 2), "\n")
  cat("\nSensitivity to the cutoff magnitude\n")
  for (cut in seq(mc, mc + 0.6, by = 0.2)) { rr <- b_mle(mags, cut); cat(sprintf("cutoff %.1f  b = %.2f +/- %.2f  n = %d\n", cut, rr["b"], rr["se"], rr["n"])) }
  cat("\nCheck on synthetic data with true b = 1.0, Mc = 4.0, 2000 events, magnitudes rounded to 0.1\n")
  sm <- round(4.0 - 0.05 + rexp(2000, rate = 1.0 * log(10)), 1)
  print(round(b_mle(sm, mc_maxcurv(sm)), 2))
  png("gr_plot.png", 800, 500); x <- seq(mc, max(mags), 0.1); N <- sapply(x, function(v) sum(mags >= v))
  plot(x, log10(N), xlab = "Magnitude", ylab = "log10 N (M >= m)"); abline(a = log10(length(mags[mags >= mc])) + r["b"] * mc, b = -r["b"], col = "red"); dev.off()
}
