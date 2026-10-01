# HG-03  Forensic DNA statistics: random match probability and the prosecutor's fallacy. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Allele frequencies below are SYNTHETIC teaching values, not a population database.
# Replace them with a published frequency table (with its terms of use checked) for real work.
demo_freqs <- list(
  L1 = c("11" = 0.15, "12" = 0.25, "13" = 0.30, "14" = 0.20, "15" = 0.10),
  L2 = c("8" = 0.10, "9" = 0.20, "10" = 0.35, "11" = 0.25, "12" = 0.10),
  L3 = c("16" = 0.20, "17" = 0.30, "18" = 0.30, "19" = 0.20),
  L4 = c("6" = 0.05, "7" = 0.30, "8" = 0.40, "9" = 0.25),
  L5 = c("20" = 0.10, "21" = 0.20, "22" = 0.30, "23" = 0.25, "24" = 0.15))

# Genotype frequency. method = "nrc410": conditional match probability of NRC II (1996) equations 4.10a and 4.10b
# (homozygote [2t+(1-t)p][3t+(1-t)p] / [(1+t)(1+2t)]; heterozygote 2[t+(1-t)pi][t+(1-t)pj] / [(1+t)(1+2t)]).
# method = "nrc44": equations 4.4a and 4.4b (homozygote p^2 + p(1-p)t; heterozygote 2 pi pj (1-t)).
# theta = 0 gives the plain Hardy-Weinberg product rule (p^2 and 2pq).
geno_freq <- function(a1, a2, freqs, theta = 0, method = c("nrc410", "nrc44")) {
  method <- match.arg(method); p1 <- freqs[[as.character(a1)]]; p2 <- freqs[[as.character(a2)]]; t <- theta
  if (method == "nrc44") return(if (a1 == a2) p1^2 + p1 * (1 - p1) * t else 2 * p1 * p2 * (1 - t))
  d <- (1 + t) * (1 + 2 * t)
  if (a1 == a2) ((2 * t + (1 - t) * p1) * (3 * t + (1 - t) * p1)) / d
  else 2 * (t + (1 - t) * p1) * (t + (1 - t) * p2) / d
}
# profile: data.frame with columns locus, a1, a2. Product rule across loci.
rmp <- function(profile, freq_list, theta = 0, method = "nrc410")
  prod(mapply(function(l, a, b) geno_freq(a, b, freq_list[[l]], theta, method), profile$locus, profile$a1, profile$a2))

# Posterior probability the suspect is the source, given RMP and prior odds 1 : (N - 1) for N possible sources.
posterior_source <- function(rmp_value, N) { prior_odds <- 1 / (N - 1); lr <- 1 / rmp_value
  odds <- prior_odds * lr; odds / (1 + odds) }
# Expected number of unrelated people in a population of size N who would match by chance.
expected_matches <- function(rmp_value, N) rmp_value * N

if (sys.nframe() == 0) {
  prof <- data.frame(locus = c("L1", "L2", "L3", "L4", "L5"), a1 = c(12, 10, 17, 8, 22), a2 = c(13, 10, 18, 9, 22))
  r0 <- rmp(prof, demo_freqs); cat(sprintf("RMP with theta = 0:    %.3e  (1 in %.0f)\n", r0, 1 / r0))
  for (th in c(0.01, 0.03)) { r <- rmp(prof, demo_freqs, th); cat(sprintf("RMP with theta = %.2f: %.3e  (1 in %.0f)\n", th, r, 1 / r)) }
  cat("\nSame profile with the NRC II 4.4 formulas, theta = 0.03:", format(rmp(prof, demo_freqs, 0.03, "nrc44"), digits = 4), "\n")
  cat("\nHeterozygote and homozygote at one locus, allele frequency 0.30 (theta 0 and 0.03)\n")
  f <- c("a" = 0.30, "b" = 0.20)
  cat(sprintf("homozygote aa: %.4f -> %.4f ; heterozygote ab: %.4f -> %.4f\n", 0.3^2, geno_freq("a", "a", as.list(f), 0.03), 2 * 0.3 * 0.2, geno_freq("a", "b", as.list(f), 0.03)))
  cat("\nProsecutor's fallacy: the RMP is not the probability the suspect is innocent\n")
  for (N in c(1e3, 1e5, 1e7)) cat(sprintf("N = %.0e possible sources: expected chance matches = %.3f, P(source | match) = %.4f\n", N, expected_matches(r0, N), posterior_source(r0, N)))
  cat("\nDatabase trawl: chance of at least one coincidental match among n unrelated people\n")
  for (n in c(1e3, 1e5, 1e6)) cat(sprintf("n = %.0e  P(>=1 match) = %.4f\n", n, 1 - (1 - r0)^n))
}
