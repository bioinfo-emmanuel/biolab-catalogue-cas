# GEN-03  Do mutations arise before selection? A Luria-Delbruck simulation. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
set.seed(1943)

# Spontaneous model: cells double each generation; each division produces a mutant with probability mu.
# Mutants then grow at the same rate as everyone else. Returns mutant counts in n_cult cultures.
spontaneous <- function(n_cult, mu, gens, n0 = 1) {
  sapply(seq_len(n_cult), function(i) { N <- n0; M <- 0
    for (g in seq_len(gens)) { new_mut <- rbinom(1, N - M, mu); M <- 2 * M + new_mut; N <- 2 * N }
    M })
}
# Induced model: no mutants until plating; each culture's count is Poisson. To compare the two models fairly the
# Poisson mean is set to a stated target, normally the mean count of the spontaneous run.
induced <- function(n_cult, target_mean) rpois(n_cult, target_mean)

vmr <- function(x) var(x) / mean(x)
# Mutation rate from the fraction of cultures with no mutants (P0 method): mu = -ln(p0) / Nf.
mu_p0 <- function(counts, Nf) { p0 <- mean(counts == 0); if (p0 == 0) NA else -log(p0) / Nf }

if (sys.nframe() == 0) {
  gens <- 20; Nf <- 2^gens; mu <- 1e-6; n <- 40
  cat(sprintf("Final cell number %.2e per culture, true mutation rate %.0e per division, %d cultures\n", Nf, mu, n))
  s <- spontaneous(n, mu, gens); i <- induced(n, mean(s))
  cat("\nSpontaneous model counts:\n"); print(sort(s)); cat(sprintf("mean %.1f  variance %.0f  variance/mean %.1f\n", mean(s), var(s), vmr(s)))
  cat("\nInduced (Poisson) model with the same mean:\n"); print(sort(i)); cat(sprintf("mean %.1f  variance %.1f  variance/mean %.2f\n", mean(i), var(i), vmr(i)))
  cat("\nVariance/mean over 200 repeat experiments of 40 cultures (induced model matched to the spontaneous mean each time)\n")
  rs <- replicate(200, { x <- spontaneous(n, mu, gens); c(vmr(x), vmr(induced(n, mean(x)))) })
  cat(sprintf("spontaneous: median %.1f (5th to 95th percentile %.1f to %.1f)\ninduced:     median %.2f (%.2f to %.2f)\n",
              median(rs[1, ]), quantile(rs[1, ], .05), quantile(rs[1, ], .95), median(rs[2, ]), quantile(rs[2, ], .05), quantile(rs[2, ], .95)))
  cat("\nMutation rate from the fraction of cultures with no mutants (100 cultures, 30 experiments)\n")
  mu2 <- 2e-7
  est <- replicate(30, mu_p0(spontaneous(100, mu2, gens), Nf)); cat(sprintf("true %.1e  mean estimate %.1e  sd %.1e\n", mu2, mean(est, na.rm = TRUE), sd(est, na.rm = TRUE)))
  cat("\nThe heavy right tail: share of the total mutant count carried by the single largest culture (300 cultures)\n")
  for (m in c(1e-7, 1e-6, 1e-5)) { x <- spontaneous(300, m, gens); cat(sprintf("mu = %.0e: cultures with no mutants %.2f; largest count %d; mean %.1f; largest culture holds %.0f%% of all mutants\n", m, mean(x == 0), max(x), mean(x), 100 * max(x) / sum(x))) }
}
