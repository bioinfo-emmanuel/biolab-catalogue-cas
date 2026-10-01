# DEV-06  French flag model of a morphogen gradient (simulation). Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
set.seed(11)

# Steady-state gradient from a source at x = 0: C(x) = C0 * exp(-x / L), L = sqrt(D / k).
gradient <- function(x, C0 = 1, L = 20) C0 * exp(-x / L)

# Two thresholds give three fates: blue above T1, white between T1 and T2, red below T2.
read_fate <- function(C, T1 = 0.6, T2 = 0.25) ifelse(C >= T1, "blue", ifelse(C >= T2, "white", "red"))

# Exact boundary positions: x = L * log(C0 / T).
boundaries <- function(C0, L, T1 = 0.6, T2 = 0.25) c(blue_white = L * log(C0 / T1), white_red = L * log(C0 / T2))

# Embryo-to-embryo variation in source strength (lognormal, coefficient of variation cv).
boundary_sd <- function(cv, L = 20, T1 = 0.6, n = 2000) {
  C0 <- exp(rnorm(n, 0, cv)); b <- L * log(C0 / T1); c(mean = mean(b), sd = sd(b))
}

if (sys.nframe() == 0) {
  x <- 0:100; L <- 20
  b <- boundaries(1, L)
  cat(sprintf("Boundaries at x = %.1f (blue/white) and x = %.1f (white/red), L = %d\n", b[1], b[2], L))
  cat("Cells per band in a 100-cell field:\n"); print(table(read_fate(gradient(x, 1, L))))
  cat("\nDoubling the source strength C0 from 1 to 2 shifts both boundaries by L * ln 2 =", round(L * log(2), 1), "\n")
  print(round(rbind(C0_1 = boundaries(1, L), C0_2 = boundaries(2, L)), 1))
  cat("\nEmbryo-to-embryo variation in C0: sd of the blue/white boundary (L = 20)\n")
  for (cv in c(0.05, 0.15, 0.30)) { r <- boundary_sd(cv); cat(sprintf("cv = %.2f  mean = %5.1f  sd = %4.2f  (approx L*cv = %.1f)\n", cv, r["mean"], r["sd"], L * cv)) }
  cat("\nWhite band width depends only on L: width = L * ln(T1/T2)\n")
  for (L2 in c(10, 20, 40)) cat(sprintf("L = %2d  width = %.1f\n", L2, diff(boundaries(1, L2))))
  cat("\nScaling: if L grows with field length (L = 0.2 * field), band proportions stay fixed\n")
  for (Fl in c(50, 100, 200)) { bb <- boundaries(1, 0.2 * Fl); cat(sprintf("field %3d  boundaries at %.0f%% and %.0f%% of the field\n", Fl, 100 * bb[1] / Fl, 100 * bb[2] / Fl)) }
}
