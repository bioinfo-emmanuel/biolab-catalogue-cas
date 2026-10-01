# GEN-02  Gene mapping from a three-point testcross (simulation). Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Genes A, B, C in that true order on the chromosome. Testcross: A B C / a b c  x  a b c / a b c.
set.seed(7)
classes <- c("ABC", "abc", "Abc", "aBC", "ABc", "abC", "AbC", "aBc")   # parental, then single-crossover classes, then double

# Class probabilities for map distances d1 (A-B) and d2 (B-C) in centimorgans and coefficient of coincidence coi.
class_probs <- function(d1, d2, coi = 1) {
  r1 <- d1 / 100; r2 <- d2 / 100; dco <- coi * r1 * r2
  sco1 <- r1 - dco; sco2 <- r2 - dco; par <- 1 - sco1 - sco2 - dco
  # Recombinant in region 1 only: aBC / Abc. Region 2 only: ABc / abC. Double: AbC / aBc.
  p <- c(ABC = par / 2, abc = par / 2, aBC = sco1 / 2, Abc = sco1 / 2, ABc = sco2 / 2, abC = sco2 / 2, AbC = dco / 2, aBc = dco / 2)
  p[classes]
}
simulate_cross <- function(n, d1, d2, coi = 1) { x <- rmultinom(1, n, class_probs(d1, d2, coi)); setNames(as.vector(x), classes) }

# Analysis from counts alone: the two most common classes are parental, the two rarest are double crossovers.
# The gene whose allele differs between a double-crossover class and its nearest parental class is the middle gene.
analyze <- function(counts) {
  o <- names(sort(counts, decreasing = TRUE)); par <- o[1:2]; dco <- o[7:8]
  chars <- function(x) strsplit(x, "")[[1]]
  ndiff <- function(x, y) sum(chars(x) != chars(y))
  p <- par[which(sapply(par, ndiff, y = dco[1]) == 1)][1]           # parental class one gene away from the double class
  middle <- which(chars(p) != chars(dco[1]))
  outer <- setdiff(1:3, middle); ord <- c(outer[1], middle, outer[2])
  n <- sum(counts)
  differs <- function(k, i) chars(k)[i] != chars(p)[i]
  recomb <- function(i, j) names(counts)[sapply(names(counts), function(k) xor(differs(k, i), differs(k, j)))]
  d1 <- sum(counts[recomb(ord[1], ord[2])]) / n; d2 <- sum(counts[recomb(ord[2], ord[3])]) / n
  dco_obs <- sum(counts[dco]) / n
  list(order = c("A", "B", "C")[ord], d_first = 100 * d1, d_second = 100 * d2, dco_obs = dco_obs,
       coi = dco_obs / (d1 * d2), interference = 1 - dco_obs / (d1 * d2))
}

if (sys.nframe() == 0) {
  cat("True map: A-B 20 cM, B-C 10 cM, coefficient of coincidence 0.5 (interference 0.5)\n")
  for (n in c(200, 1000, 10000)) {
    cnt <- simulate_cross(n, 20, 10, 0.5); a <- analyze(cnt)
    cat(sprintf("\nn = %d\n", n)); print(cnt)
    cat(sprintf("order %s | first interval %.1f cM | second interval %.1f cM | coincidence %.2f | interference %.2f\n", paste(a$order, collapse = "-"), a$d_first, a$d_second, a$coi, a$interference))
  }
  cat("\nSpread of estimates over 500 experiments of 500 progeny\n")
  est <- t(replicate(500, { a <- analyze(simulate_cross(500, 20, 10, 0.5)); c(a$d_first, a$d_second) }))
  print(round(rbind(mean = colMeans(est), sd = apply(est, 2, sd)), 2))
  cat("\nMap distance versus recombinant fraction: intervals of 40 and 20 cM, no interference\n")
  p <- class_probs(40, 20, 1)
  cat("Map distance A-C = 60 cM, but the observed A-C recombinant fraction is", round(sum(p[c("aBC", "Abc", "ABc", "abC")]), 3),
      "because double crossovers restore the parental combination of A and C.\n")
}
