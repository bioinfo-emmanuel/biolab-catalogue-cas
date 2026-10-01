# HG-05  Heritability from twins: simulating the ACE model. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Phenotype P = a*A + c*C + e*E, with A additive genetic, C shared environment, E unique environment.
# MZ twins share all of A, DZ twins share half; both share C.
set.seed(101)
sim_twins <- function(n_pairs, h2, c2, zyg = c("MZ", "DZ"), extra_c_mz = 0) {
  zyg <- match.arg(zyg)
  if (zyg == "MZ") c2 <- c2 + extra_c_mz   # extra shared environment for MZ only (equal-environment violation)
  e2 <- 1 - h2 - c2; r_A <- if (zyg == "MZ") 1 else 0.5
  A_shared <- rnorm(n_pairs); A1 <- sqrt(r_A) * A_shared + sqrt(1 - r_A) * rnorm(n_pairs)
  A2 <- sqrt(r_A) * A_shared + sqrt(1 - r_A) * rnorm(n_pairs)
  C <- rnorm(n_pairs)
  data.frame(zyg = zyg, p1 = sqrt(h2) * A1 + sqrt(c2) * C + sqrt(e2) * rnorm(n_pairs),
                        p2 = sqrt(h2) * A2 + sqrt(c2) * C + sqrt(e2) * rnorm(n_pairs))
}
# Falconer estimates from twin correlations.
falconer <- function(d) {
  rMZ <- with(d[d$zyg == "MZ", ], cor(p1, p2)); rDZ <- with(d[d$zyg == "DZ", ], cor(p1, p2))
  c(rMZ = rMZ, rDZ = rDZ, h2 = 2 * (rMZ - rDZ), c2 = 2 * rDZ - rMZ, e2 = 1 - rMZ)
}
# Bootstrap standard error by resampling pairs within zygosity.
boot_se <- function(d, B = 500) {
  est <- replicate(B, { s <- do.call(rbind, lapply(split(d, d$zyg), function(g) g[sample(nrow(g), replace = TRUE), ])); falconer(s) })
  apply(est, 1, sd)
}
if (sys.nframe() == 0) {
  cat("True h2 = 0.5, c2 = 0.2, e2 = 0.3\n")
  for (n in c(50, 200, 1000)) {
    d <- rbind(sim_twins(n, 0.5, 0.2, "MZ"), sim_twins(n, 0.5, 0.2, "DZ"))
    cat(sprintf("\n%d pairs of each type\n", n)); print(round(rbind(estimate = falconer(d), se = boot_se(d, 300)), 2))
  }
  cat("\nRepeat 300 studies with 100 pairs per type: spread of the h2 estimate\n")
  h <- replicate(300, falconer(rbind(sim_twins(100, 0.5, 0.2, "MZ"), sim_twins(100, 0.5, 0.2, "DZ")))["h2"])
  cat(sprintf("mean = %.2f  sd = %.2f  range = %.2f to %.2f\n", mean(h), sd(h), min(h), max(h)))
  cat("\nEqual-environment assumption broken: true h2 = 0.3, c2 = 0.1, MZ pairs get 0.1 extra shared environment (5000 pairs each)\n")
  d <- rbind(sim_twins(5000, 0.3, 0.1, "MZ", extra_c_mz = 0.1), sim_twins(5000, 0.3, 0.1, "DZ"))
  print(round(falconer(d), 2)); cat("Falconer h2 should be near 0.3 but is pushed up by about 2 * 0.1\n")
}
