# DEV-01  Temperature and developmental rate in fruit flies: degree-days. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript dev01_degree_days.R [vials.csv]
# vials.csv: vial, temp_c (measured mean temperature of the vial), egg_to_adult_days (mean or median day of adult emergence from the day eggs were laid), setting_c (optional: the nominal temperature group)
# Model: rate = 1 / days = (T - T0) / K, so the rate is a straight line in T; T0 is the lower threshold and K the thermal constant in degree-days.
set.seed(1)
fit_dd <- function(d) { d$rate <- 1 / d$egg_to_adult_days; f <- lm(rate ~ temp_c, d); a <- coef(f)[[1]]; b <- coef(f)[[2]]; list(fit = f, T0 = -a / b, K = 1 / b) }
boot_dd <- function(d, B = 2000) { r <- t(replicate(B, { x <- d[sample(nrow(d), replace = TRUE), ]; if (length(unique(x$temp_c)) < 2) c(NA, NA) else { f <- fit_dd(x); c(f$T0, f$K) } })); r <- r[complete.cases(r), ]; apply(r, 2, quantile, c(0.025, 0.975)) }
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) d <- read.csv(args[1]) else {
    T0 <- 9; K <- 170; temps <- c(19.5, 22.5, 26.0, 29.0); d <- do.call(rbind, lapply(temps, function(t) data.frame(vial = paste0("v", t, "_", 1:4), setting_c = t, temp_c = t + rnorm(4, 0, 0.3), egg_to_adult_days = round(K / (t - T0) * rnorm(4, 1, 0.04), 1))))
    cat("Practice data: 4 temperatures x 4 vials with invented values (true lower threshold 9 C, thermal constant 170 degree-days)\n\n") }
  grp <- if ("setting_c" %in% names(d)) d$setting_c else round(d$temp_c); m <- do.call(rbind, lapply(split(d, grp), function(x) data.frame(temp_c = mean(x$temp_c), n_vials = nrow(x), mean_days = mean(x$egg_to_adult_days), sd_days = sd(x$egg_to_adult_days), mean_rate = mean(1 / x$egg_to_adult_days)))); print(format(m, digits = 3), row.names = FALSE)
  f <- fit_dd(d); cat(sprintf("\nLinear model of rate on temperature: slope %.5f per C per day, R squared %.3f\nLower threshold T0 = %.1f C, thermal constant K = %.0f degree-days\n", coef(f$fit)[2], summary(f$fit)$r.squared, f$T0, f$K))
  b <- boot_dd(d); cat(sprintf("Bootstrap 95%% intervals: T0 %.1f to %.1f C; K %.0f to %.0f degree-days\n", b[1, 1], b[2, 1], b[1, 2], b[2, 2]))
  cat("\nPrediction: development time at 25 C =", round(f$K / (25 - f$T0), 1), "days\n")
  cat("Extrapolation warning: T0 lies well below the tested temperatures, so it is an extrapolation and is sensitive to small errors in the rate at the coldest temperature.\n")
  cat("Check: K = 170 and T0 = 9 give", round(170 / (27 - 9), 1), "days at 27 C, and 170 / (18 - 9) =", round(170 / 9, 1), "days at 18 C\n")
  cat("\nEffect of the number of temperatures: 3 versus 2 (removing the coldest):\n"); f2 <- fit_dd(d[d$temp_c > min(d$temp_c) + 1, ]); cat(sprintf("T0 %.1f, K %.0f from the warmer three; from all four T0 %.1f, K %.0f\n", f2$T0, f2$K, f$T0, f$K))
}
