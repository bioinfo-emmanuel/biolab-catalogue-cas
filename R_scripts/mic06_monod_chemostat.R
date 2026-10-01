# MIC-06  Monod growth and the chemostat (simulation). Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
monod <- function(S, mu_max, Ks) mu_max * S / (Ks + S)

# Chemostat with dilution rate D, feed concentration S0, yield Y.
# dX/dt = (mu - D) X ; dS/dt = D (S0 - S) - mu X / Y   (solved with RK4)
chemostat <- function(D, mu_max = 0.6, Ks = 0.2, Y = 0.5, S0 = 5, X0 = 0.05, S_init = 5, tmax = 100, dt = 0.05) {
  f <- function(u) { mu <- monod(u[2], mu_max, Ks)
    c((mu - D) * u[1], D * (S0 - u[2]) - mu * u[1] / Y) }
  n <- round(tmax / dt); out <- matrix(NA, n + 1, 3); u <- c(X0, S_init); out[1, ] <- c(0, u)
  for (i in 1:n) { k1 <- f(u); k2 <- f(u + dt / 2 * k1); k3 <- f(u + dt / 2 * k2); k4 <- f(u + dt * k3)
    u <- u + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4); out[i + 1, ] <- c(i * dt, u) }
  colnames(out) <- c("t", "X", "S"); as.data.frame(out)
}

steady_state <- function(D, mu_max = 0.6, Ks = 0.2, Y = 0.5, S0 = 5) {
  if (D >= mu_max * S0 / (Ks + S0)) return(c(S = S0, X = 0))   # washout
  S <- D * Ks / (mu_max - D); c(S = S, X = Y * (S0 - S))
}

# Fit mu_max and Ks from batch-culture specific growth rates by nonlinear least squares.
fit_monod <- function(S, mu) nls(mu ~ mu_max * S / (Ks + S), start = list(mu_max = max(mu) * 1.2, Ks = median(S)))

if (sys.nframe() == 0) {
  cat("Steady state vs dilution rate D (mu_max 0.6, Ks 0.2, Y 0.5, S0 5)\n")
  Ds <- c(0.1, 0.3, 0.5, 0.55, 0.6)
  print(round(cbind(D = Ds, t(sapply(Ds, steady_state))), 3))
  cat(sprintf("\nCritical dilution rate (washout above): %.3f per h\n", 0.6 * 5 / (0.2 + 5)))
  d <- chemostat(0.3); cat("Simulated X, S at t = 100 with D = 0.3:", round(as.numeric(tail(d, 1)[2:3]), 3), "\n")
  d2 <- chemostat(0.58); cat("Simulated X at t = 100 with D = 0.58 (above critical):", round(tail(d2$X, 1), 5), "\n")
  set.seed(5); S <- c(0.05, 0.1, 0.2, 0.5, 1, 2, 5); mu <- monod(S, 0.6, 0.2) * exp(rnorm(7, 0, 0.04))
  fit <- fit_monod(S, mu); cat("\nFit of noisy batch data:\n"); print(round(coef(summary(fit))[, 1:2], 3))
}
