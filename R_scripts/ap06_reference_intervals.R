# AP-06  Reference intervals and normal variation from a population health dataset. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript ap06_reference_intervals.R [data.csv]
# data.csv columns: age (years), sex (F or M), value (one measurement, for example systolic blood pressure).
# With no argument the script simulates 1500 adults so it runs offline. Use only de-identified open data.
set.seed(2006)

simulate_health <- function(n = 1500) { age <- round(runif(n, 18, 80)); sex <- sample(c("F", "M"), n, TRUE)
  sbp <- 98 + 0.55 * age + ifelse(sex == "M", 6, 0) + rnorm(n, 0, 11); data.frame(age = age, sex = sex, value = round(sbp)) }
ref_int <- function(x, level = 0.95) { a <- (1 - level) / 2; quantile(x, c(a, 1 - a), names = FALSE, type = 7) }
boot_limits <- function(x, B = 1000) { b <- replicate(B, ref_int(sample(x, replace = TRUE))); rbind(lower_limit_ci = quantile(b[1, ], c(.05, .95)), upper_limit_ci = quantile(b[2, ], c(.05, .95))) }

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE); d <- if (length(args) >= 1) read.csv(args[1]) else { cat("Practice data: 1500 simulated adults, aged 18 to 80\n"); simulate_health() }
  d$age_group <- cut(d$age, c(17, 34, 49, 64, 80), labels = c("18-34", "35-49", "50-64", "65-80"))
  cat(sprintf("\nAll adults (n = %d): mean %.1f, sd %.1f; 95%% reference interval %.0f to %.0f\n", nrow(d), mean(d$value), sd(d$value), ref_int(d$value)[1], ref_int(d$value)[2]))
  cat("\nReference intervals (2.5th to 97.5th percentile) by age group and sex\n")
  tab <- aggregate(value ~ age_group + sex, d, function(x) c(n = length(x), lower = ref_int(x)[1], median = median(x), upper = ref_int(x)[2])); tab <- do.call(data.frame, tab); names(tab) <- c("age_group", "sex", "n", "lower", "median", "upper"); print(tab, digits = 4, row.names = FALSE)
  fit <- lm(value ~ age + sex, d); cat("\nLinear model of the measurement on age and sex:\n"); print(round(summary(fit)$coefficients, 3))
  cat(sprintf("Each extra year of age changes the measurement by %.2f units (95%% CI %.2f to %.2f). Residual sd %.1f; R squared %.2f.\n", coef(fit)["age"], confint(fit)["age", 1], confint(fit)["age", 2], sigma(fit), summary(fit)$r.squared))
  cat("\nHow sure are the limits? Bootstrap 90% intervals for the limits of the whole sample, by sample size\n")
  for (n in c(40, 120, 400)) { x <- sample(d$value, n); b <- boot_limits(x, 500); cat(sprintf("n = %3d  lower limit %.0f (%.0f to %.0f)   upper limit %.0f (%.0f to %.0f)\n", n, ref_int(x)[1], b[1, 1], b[1, 2], ref_int(x)[2], b[2, 1], b[2, 2])) }
  cat("\nA reference interval flags 5% of healthy people by design.\n")
  cat("Chance that a healthy person has at least one 'abnormal' result among 10 independent tests: expected", round(1 - 0.95^10, 3), ", simulated", round(mean(replicate(20000, any(runif(10) < 0.05))), 3), "\n")
  cat("\nUsing one interval for everyone versus age-specific intervals: share of the oldest group outside the whole-sample interval\n")
  ri <- ref_int(d$value); old <- d$value[d$age >= 65]; cat(sprintf("%.1f%% of adults aged 65 to 80 fall outside the whole-sample interval (5%% expected if age did not matter)\n", 100 * mean(old < ri[1] | old > ri[2])))
}
