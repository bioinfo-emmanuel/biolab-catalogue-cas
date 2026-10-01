# DEV-04  Avian embryo staging from eggs of known incubation age. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript dev04_embryo_staging.R [embryos.csv]
# embryos.csv: egg, hours_incubated, observer, somites (count of somite pairs), hh_stage (Hamburger-Hamilton stage you assigned)
# Reference timings (Hamburger and Hamilton 1951; table as reproduced by UNSW Embryology): stage 4 at 18-19 h, 8 at 26-29 h (4 somites), 10 at 33-38 h (10 somites),
# 12 at 45-49 h (16 somites), 14 at 50-53 h (22 somites), 16 at 51-56 h (26-28 somites), 18 at 3 days (30-36 somites; 72 h used here). Eggs incubated at 37.5 C.
set.seed(4)
ref <- data.frame(stage = c(4, 8, 10, 12, 14, 16, 18), hours_lo = c(18, 26, 33, 45, 50, 51, 72), hours_hi = c(19, 29, 38, 49, 53, 56, 72), somites = c(0, 4, 10, 16, 22, 27, 33))
expected_stage <- function(h) approx(rowMeans(ref[, c("hours_lo", "hours_hi")]), ref$stage, xout = h, rule = 2)$y
cohen_kappa <- function(a, b) { lv <- sort(union(a, b)); tab <- table(factor(a, lv), factor(b, lv)); n <- sum(tab); po <- sum(diag(tab)) / n; pe <- sum(rowSums(tab) * colSums(tab)) / n^2; if (pe == 1) 1 else (po - pe) / (1 - pe) }
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) d <- read.csv(args[1]) else {
    hrs <- rep(c(24, 36, 48, 60, 72), each = 4); n <- length(hrs); som <- pmax(0, round(rnorm(n, approx(rowMeans(ref[, 2:3]), ref$somites, xout = hrs, rule = 2)$y, 2)))
    st <- function(s, h) ref$stage[which.min(abs(ref$somites - s))]; a <- data.frame(egg = 1:n, hours_incubated = hrs, observer = "A", somites = som, hh_stage = sapply(seq_len(n), function(i) st(som[i], hrs[i])))
    b <- a; b$observer <- "B"; b$somites <- pmax(0, som + sample(-1:1, n, TRUE)); b$hh_stage <- sapply(seq_len(n), function(i) st(b$somites[i], hrs[i])); d <- rbind(a, b)
    cat("Practice data: 20 invented eggs opened after 24 to 72 hours, staged by two observers from somite counts; the counts are invented\n\n") }
  a <- d[d$observer == unique(d$observer)[1], ]
  cat("Mean somite number and stage by hours of incubation:\n"); print(round(aggregate(cbind(somites, hh_stage) ~ hours_incubated, a, mean), 1), row.names = FALSE)
  f <- lm(somites ~ hours_incubated, a); cat(sprintf("\nSomite pairs formed per hour of incubation (regression): %.2f (95%% CI %.2f to %.2f); about %.1f hours per somite\n", coef(f)[2], confint(f)[2, 1], confint(f)[2, 2], 1 / coef(f)[2]))
  a$expected <- expected_stage(a$hours_incubated); cat(sprintf("Mean difference between the stage you assigned and the stage expected from the incubation time: %+.1f stages (sd %.1f)\n", mean(a$hh_stage - a$expected), sd(a$hh_stage - a$expected)))
  cat(sprintf("Within-time spread: at 48 hours the stages assigned ranged from %d to %d, so incubation time predicts stage only approximately\n", min(a$hh_stage[a$hours_incubated == 48]), max(a$hh_stage[a$hours_incubated == 48])))
  w <- reshape(d[, c("egg", "observer", "hh_stage")], idvar = "egg", timevar = "observer", direction = "wide"); cat(sprintf("\nObserver agreement on stage: exact agreement %.0f%%, within two stage numbers %.0f%%, kappa %.2f\n", 100 * mean(w[[2]] == w[[3]]), 100 * mean(abs(w[[2]] - w[[3]]) <= 2), cohen_kappa(w[[2]], w[[3]])))
  cat("Stages here are the reference stages 4, 8, 10, 12, 14, 16 and 18 only; real staging uses every stage from a photograph series.\n")
  cat("\nReference table used in this script (hours are ranges from the published table):\n"); print(ref, row.names = FALSE)
}
