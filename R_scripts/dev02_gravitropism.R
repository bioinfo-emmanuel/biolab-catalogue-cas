# DEV-02  Gravitropism in seedling roots: kinetics of curvature. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript dev02_gravitropism.R [angles.csv]
# angles.csv: treatment, seedling, hours, angle_deg   (angle of the growing root tip from its original direction after the dish was turned through 90 degrees; 0 = no bending)
set.seed(2)
lag_and_rate <- function(x) { x <- x[order(x$hours), ]; if (nrow(x) < 4) return(c(lag_h = NA, max_rate = NA, final = NA)); s <- diff(x$angle_deg) / diff(x$hours); i <- which.max(s)
  lag <- if (any(x$angle_deg >= 10)) approx(x$angle_deg, x$hours, xout = 10, ties = "ordered")$y else NA; c(lag_h = lag, max_rate = s[i], final = tail(x$angle_deg, 1)) }
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) d <- read.csv(args[1]) else {
    hrs <- c(0, 1, 2, 3, 4, 6, 8, 12, 24); spec <- list(room = c(lag = 1.2, tmax = 4.0, max = 86), cool = c(lag = 2.2, tmax = 8.0, max = 78))
    d <- do.call(rbind, lapply(names(spec), function(tr) do.call(rbind, lapply(1:10, function(s) { p <- spec[[tr]]; a <- p[["max"]] / (1 + exp(-(hrs - p[["lag"]] - p[["tmax"]] / 2 * rnorm(1, 1, 0.1)) / (p[["tmax"]] / 5))); a[1] <- 0; data.frame(treatment = tr, seedling = paste0(tr, s), hours = hrs, angle_deg = round(pmin(90, pmax(0, a + rnorm(length(hrs), 0, 3))))) }))))
    cat("Practice data: 10 invented seedlings in each of two invented conditions (room temperature and cooler), angles at 9 times over 24 hours\n\n") }
  m <- aggregate(angle_deg ~ treatment + hours, d, function(x) c(mean = mean(x), sd = sd(x))); m <- do.call(data.frame, m); m[, 3:4] <- round(m[, 3:4], 1); cat("Mean angle by time:\n"); print(reshape(m[, c("treatment", "hours", "angle_deg.mean")], idvar = "hours", timevar = "treatment", direction = "wide"), row.names = FALSE)
  k <- do.call(rbind, lapply(split(d, d$seedling), function(x) data.frame(treatment = x$treatment[1], seedling = x$seedling[1], t(lag_and_rate(x))))); res <- aggregate(cbind(lag_h, max_rate, final) ~ treatment, k, function(v) c(mean = mean(v), sd = sd(v))); res <- do.call(data.frame, res); res[, -1] <- round(res[, -1], 2); cat("\nPer-seedling summaries (lag to 10 degrees, steepest bending rate in degrees per hour, final angle):\n"); print(res, row.names = FALSE)
  for (v in c("lag_h", "max_rate", "final")) { w <- wilcox.test(k[[v]] ~ k$treatment, exact = FALSE); cat(sprintf("%s: Wilcoxon p = %.4f\n", v, w$p.value)) }
  cat("\nCheck: the angle of a root that grows 4 mm horizontally and 4 mm downward after the turn is", round(atan2(4, 4) * 180 / pi), "degrees from its original direction\n")
  cat("The seedling is the unit of replication; readings on the same seedling at different times are not independent.\n")
}
