# ENV-05  Radioactive decay and dating uncertainty (simulation). Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: source("env05_decay_dating.R")  or  Rscript env05_decay_dating.R
set.seed(2026)

# 1. Decay of N0 atoms: each atom decays with probability p = 1 - exp(-lambda * dt) per step.
simulate_decay <- function(N0, half_life, dt, steps) {
  lambda <- log(2) / half_life
  p <- 1 - exp(-lambda * dt)
  N <- numeric(steps + 1); N[1] <- N0
  for (i in 1:steps) N[i + 1] <- N[i] - rbinom(1, N[i], p)
  data.frame(t = (0:steps) * dt, N = N)
}

# 2. Estimate half-life from one run by regressing ln(N) on time (N > 0 only).
estimate_half_life <- function(d) {
  d <- d[d$N > 0, ]
  fit <- lm(log(N) ~ t, data = d)
  -log(2) / coef(fit)[["t"]]
}

# 3. Age from a measured parent fraction f = N/N0: t = -half_life * log2(f).
age_from_fraction <- function(f, half_life) -half_life * log2(f)

# 4. Uncertainty of the age when the count is Poisson: sd(N) = sqrt(N).
age_ci <- function(N_obs, N0, half_life, level = 0.95) {
  z <- qnorm(1 - (1 - level) / 2)
  f <- N_obs / N0
  f_lo <- max(N_obs - z * sqrt(N_obs), 1) / N0
  f_hi <- (N_obs + z * sqrt(N_obs)) / N0
  c(age = age_from_fraction(f, half_life),
    lower = age_from_fraction(min(f_hi, 1), half_life),
    upper = age_from_fraction(f_lo, half_life))
}

if (sys.nframe() == 0) {
  half_life <- 100                      # arbitrary time units
  cat("Effect of starting size N0 on the half-life estimate (200 runs each)\n")
  for (N0 in c(50, 500, 5000)) {
    est <- replicate(200, estimate_half_life(simulate_decay(N0, half_life, 10, 40)))
    cat(sprintf("N0 = %5d  mean = %6.1f  sd = %5.1f  cv = %4.1f%%\n",
                N0, mean(est), sd(est), 100 * sd(est) / mean(est)))
  }
  cat("\nAge and 95% interval from a sample with 250 of 1000 parent atoms left\n")
  print(round(age_ci(250, 1000, half_life), 1))
  cat("\nSame fraction (0.25) with 25 of 100 and 2500 of 10000 atoms\n")
  print(round(age_ci(25, 100, half_life), 1)); print(round(age_ci(2500, 10000, half_life), 1))
  d <- simulate_decay(1000, half_life, 10, 40)
  png("decay_curve.png", 800, 500); plot(d$t, d$N, type = "b", xlab = "Time", ylab = "Atoms remaining")
  curve(1000 * 2^(-x / half_life), add = TRUE, col = "red"); dev.off()
  cat("\nSaved decay_curve.png\n")
}
