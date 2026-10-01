# EVO-05  Testing Hardy-Weinberg proportions, and the Wahlund effect. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript evo05_hardy_weinberg.R [counts.csv]
# counts.csv columns: population, AA, Aa, aa   (genotype counts at one biallelic locus).
# With no argument the script simulates practice populations so it runs offline.
set.seed(315)

hw_test <- function(n_AA, n_Aa, n_aa) {
  n <- n_AA + n_Aa + n_aa; p <- (2 * n_AA + n_Aa) / (2 * n); q <- 1 - p
  exp <- n * c(p^2, 2 * p * q, q^2); obs <- c(n_AA, n_Aa, n_aa)
  chi <- sum((obs - exp)^2 / exp)
  c(n = n, p = p, He = 2 * p * q, Ho = n_Aa / n, F = 1 - (n_Aa / n) / (2 * p * q), chisq = chi, p_value = pchisq(chi, 1, lower.tail = FALSE))
}
# Exact test by enumerating heterozygote counts given the allele counts (Wigginton et al. 2005 idea, done by brute force).
hw_exact <- function(n_AA, n_Aa, n_aa) {
  n <- n_AA + n_Aa + n_aa; nA <- 2 * n_AA + n_Aa; na <- 2 * n_aa + n_Aa
  hs <- seq(nA %% 2, min(nA, na), by = 2)
  lp <- sapply(hs, function(h) { nAA <- (nA - h) / 2; nAa <- (na - h) / 2
    lfactorial(n) - lfactorial(nAA) - lfactorial(h) - lfactorial(nAa) + h * log(2) - (lfactorial(2 * n) - lfactorial(nA) - lfactorial(na)) })
  pr <- exp(lp - max(lp)); pr <- pr / sum(pr); pobs <- pr[hs == n_Aa]
  sum(pr[pr <= pobs + 1e-12])
}
simulate_pop <- function(n, p, F = 0) { q <- 1 - p
  pr <- c(p^2 + F * p * q, 2 * p * q * (1 - F), q^2 + F * p * q); x <- rmultinom(1, n, pr); setNames(as.vector(x), c("AA", "Aa", "aa")) }

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) dat <- read.csv(args[1]) else {
    dat <- data.frame(population = c("A", "B", "C", "D"), rbind(simulate_pop(200, 0.8, 0), simulate_pop(200, 0.2, 0), simulate_pop(200, 0.5, 0), simulate_pop(200, 0.5, 0.25)))
    cat("Practice data: four populations of 200; the fourth has a true inbreeding coefficient F = 0.25\n") }
  res <- t(sapply(seq_len(nrow(dat)), function(i) hw_test(dat$AA[i], dat$Aa[i], dat$aa[i])))
  res <- data.frame(population = dat$population, round(res, 3), exact_p = sapply(seq_len(nrow(dat)), function(i) round(hw_exact(dat$AA[i], dat$Aa[i], dat$aa[i]), 3)))
  print(res, row.names = FALSE)
  cat("\nWahlund effect: pool populations A and B (frequencies 0.8 and 0.2, each in Hardy-Weinberg proportions)\n")
  pooled <- unname(colSums(dat[dat$population %in% c("A", "B"), c("AA", "Aa", "aa")])); pr <- hw_test(pooled[1], pooled[2], pooled[3])
  print(round(pr, 3))
  cat("\nPower: chance of rejecting HW at 0.05 when F = 0.25, by sample size (500 simulations each)\n")
  for (n in c(20, 50, 100, 200, 500)) { rej <- mean(replicate(500, { x <- simulate_pop(n, 0.5, 0.25); hw_test(x[[1]], x[[2]], x[[3]])["p_value"] < 0.05 })); cat(sprintf("n = %3d  power = %.2f\n", n, rej)) }
  cat("\nFalse positive rate when HW holds (n = 100, 2000 simulations): ")
  cat(round(mean(replicate(2000, { x <- simulate_pop(100, 0.3, 0); hw_test(x[[1]], x[[2]], x[[3]])["p_value"] < 0.05 })), 3), "\n")
}
