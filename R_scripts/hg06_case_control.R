# HG-06  Case-control association and population stratification. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
set.seed(2027)
# Two subpopulations differ in allele frequency (Fst-like difference) and in baseline disease risk.
sim_study <- function(n_cases = 500, n_controls = 500, p1 = 0.2, p2 = 0.5, risk1 = 0.02, risk2 = 0.06,
                      frac_pop2_cases = NULL, causal_or = 1) {
  # Case-control sampling: draw cases and controls from a population where pop2 has higher risk.
  pop_frac2 <- 0.5
  w_case <- c(1 - pop_frac2, pop_frac2) * c(risk1, risk2); w_case <- w_case / sum(w_case)
  w_ctrl <- c(1 - pop_frac2, pop_frac2) * (1 - c(risk1, risk2)); w_ctrl <- w_ctrl / sum(w_ctrl)
  pop_case <- rbinom(n_cases, 1, w_case[2]) + 1; pop_ctrl <- rbinom(n_controls, 1, w_ctrl[2]) + 1
  pop <- c(pop_case, pop_ctrl); status <- rep(c(1, 0), c(n_cases, n_controls)); p <- c(p1, p2)[pop]
  g <- rbinom(length(pop), 2, p)   # allele count with no causal effect
  data.frame(status = status, pop = pop, g = g)
}
allelic_test <- function(d) {
  a1 <- tapply(d$g, d$status, sum); n <- tapply(d$g, d$status, length) * 2
  m <- rbind(cases = c(a1["1"], n["1"] - a1["1"]), controls = c(a1["0"], n["0"] - a1["0"]))
  ct <- chisq.test(m, correct = FALSE)
  c(OR = (m[1, 1] * m[2, 2]) / (m[1, 2] * m[2, 1]), chisq = unname(ct$statistic), p = ct$p.value)
}
# Stratified analysis: Cochran-Mantel-Haenszel over populations (allele-count 2x2 tables).
cmh_test <- function(d) {
  tabs <- lapply(1:2, function(k) { s <- d[d$pop == k, ]
    a1 <- tapply(s$g, s$status, sum); n <- tapply(s$g, s$status, length) * 2
    matrix(c(a1["1"], n["1"] - a1["1"], a1["0"], n["0"] - a1["0"]), 2, byrow = TRUE) })
  num <- sum(sapply(tabs, function(t) t[1, 1] * t[2, 2] / sum(t))); den <- sum(sapply(tabs, function(t) t[1, 2] * t[2, 1] / sum(t)))
  arr <- array(unlist(tabs), c(2, 2, 2)); mt <- mantelhaen.test(arr, correct = FALSE)
  c(OR_MH = num / den, chisq = unname(mt$statistic), p = mt$p.value)
}
# Genomic inflation factor from many null tests: median(chisq) / 0.456.
lambda_gc <- function(chisq) median(chisq) / qchisq(0.5, 1)

if (sys.nframe() == 0) {
  d <- sim_study()
  cat("One study, no causal effect, populations differ in allele frequency and risk\n")
  print(round(rbind(naive = allelic_test(d), stratified_CMH = cmh_test(d)), 4))
  cat("\n1000 null markers: false positive rate at p < 0.05 and inflation factor\n")
  res <- t(replicate(1000, { dd <- sim_study(); c(naive = allelic_test(dd)["chisq"], cmh = cmh_test(dd)["chisq"]) }))
  colnames(res) <- c("naive", "cmh")
  fp <- function(x) mean(1 - pchisq(x, 1) < 0.05)
  cat(sprintf("naive: FPR = %.2f, lambda = %.2f\ncmh:   FPR = %.2f, lambda = %.2f\n", fp(res[, 1]), lambda_gc(res[, 1]), fp(res[, 2]), lambda_gc(res[, 2])))
  cat("\nNo stratification (same allele frequency and risk in both populations):\n")
  res0 <- replicate(1000, allelic_test(sim_study(p1 = 0.35, p2 = 0.35, risk1 = 0.04, risk2 = 0.04))["chisq"])
  cat(sprintf("naive: FPR = %.2f, lambda = %.2f\n", fp(res0), lambda_gc(res0)))
}
