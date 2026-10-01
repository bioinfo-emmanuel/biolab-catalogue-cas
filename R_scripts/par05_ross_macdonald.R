# PAR-05  Ross-Macdonald model of vector-borne transmission (simulation). Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# State: x = fraction of humans infected, y = fraction of mosquitoes infectious.
# dx/dt = m a b y (1 - x) - r x
# dy/dt = a c x (1 - y) exp(-g n) - g y        (Macdonald form; n = extrinsic incubation period)
# Parameter values are illustrative and must be replaced from a cited source for any real use.
default_pars <- list(m = 4, a = 0.25, b = 0.4, c = 0.4, g = 0.15, n = 10, r = 0.01)

R0 <- function(p) with(p, m * a^2 * b * c * exp(-g * n) / (r * g))

rk4 <- function(f, u, dt, tmax) {
  n <- round(tmax / dt); out <- matrix(NA, n + 1, length(u) + 1); out[1, ] <- c(0, u)
  for (i in 1:n) { k1 <- f(u); k2 <- f(u + dt / 2 * k1); k3 <- f(u + dt / 2 * k2); k4 <- f(u + dt * k3)
    u <- u + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4); out[i + 1, ] <- c(i * dt, u) }
  out
}
simulate <- function(p, x0 = 0.01, y0 = 0, tmax = 1500, dt = 0.5) {
  f <- function(u) with(p, c(m * a * b * u[2] * (1 - u[1]) - r * u[1],
                             a * c * u[1] * (1 - u[2]) * exp(-g * n) - g * u[2]))
  o <- as.data.frame(rk4(f, c(x0, y0), dt, tmax)); names(o) <- c("t", "x", "y"); o
}
# Effect of vector control: change one parameter and recompute R0 and the endemic level.
scenario <- function(p, ...) { q <- modifyList(p, list(...)); s <- simulate(q)
  c(R0 = R0(q), x_end = tail(s$x, 1), y_end = tail(s$y, 1)) }

if (sys.nframe() == 0) {
  p <- default_pars
  cat(sprintf("Baseline R0 = %.2f\n", R0(p)))
  res <- rbind(baseline = scenario(p), bednets_a_x0.5 = scenario(p, a = 0.125),
               adulticide_g_x2 = scenario(p, g = 0.30), larval_control_m_x0.25 = scenario(p, m = 1))
  print(round(res, 3))
  cat("\nThreshold: mosquito density m at which R0 = 1\n")
  m_crit <- uniroot(function(m) R0(modifyList(p, list(m = m))) - 1, c(0.01, 100))$root
  cat(sprintf("m_crit = %.2f mosquitoes per person\n", m_crit))
  cat("\nSensitivity of R0 to a 10% rise in each parameter\n")
  base <- R0(p); sens <- sapply(names(p), function(k) { q <- p; q[[k]] <- q[[k]] * 1.1; R0(q) / base - 1 })
  print(round(100 * sens, 1))
}
