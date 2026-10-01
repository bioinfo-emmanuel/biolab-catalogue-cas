# ZOO-03  Choice experiments on invertebrate behaviour with proper replication. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript zoo03_choice_test.R [choices.csv]
# choices.csv: trial, animal, arena (choice or control), start_side (A or B), side_chosen (A or B), seconds_on_A, seconds_on_B
# Condition A is the treatment side (for example damp). Control arenas have no difference between sides and test for side bias.
set.seed(3)
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) d <- read.csv(args[1]) else {
    n <- 30; mk <- function(arena, pA, tA) { ch <- ifelse(runif(n) < pA, "A", "B"); on_a <- ifelse(ch == "A", runif(n, tA, 300), runif(n, 0, 300 - tA)); data.frame(trial = 1:n, animal = paste0(arena, "_", 1:n), arena = arena, start_side = sample(c("A", "B"), n, TRUE), side_chosen = ch, seconds_on_A = round(on_a), seconds_on_B = round(300 - on_a)) }
    d <- rbind(mk("choice", 0.80, 120), mk("control", 0.50, 120)); cat("Practice data: 30 woodlice-like animals in each arena, each used once; 300 s trials. The 80 percent preference is invented.\n\n") }
  ch <- d[d$arena == "choice", ]; ct <- d[d$arena == "control", ]
  b <- binom.test(sum(ch$side_chosen == "A"), nrow(ch), 0.5); cat(sprintf("Choice arena: %d of %d animals chose A (%.0f%%), 95%% CI %.2f to %.2f, exact binomial p = %.4f\n", sum(ch$side_chosen == "A"), nrow(ch), 100 * mean(ch$side_chosen == "A"), b$conf.int[1], b$conf.int[2], b$p.value))
  bc <- binom.test(sum(ct$side_chosen == "A"), nrow(ct), 0.5); cat(sprintf("Control arena (no gradient): %d of %d chose side A, p = %.3f (a significant result here would point to a side bias in the set-up)\n", sum(ct$side_chosen == "A"), nrow(ct), bc$p.value))
  f <- fisher.test(table(d$arena, d$side_chosen)); cat(sprintf("Choice versus control: Fisher exact p = %.4f, odds ratio %.1f\n", f$p.value, f$estimate))
  w <- wilcox.test(ch$seconds_on_A, ct$seconds_on_A, exact = FALSE); cat(sprintf("Time on side A, choice versus control: medians %.0f and %.0f s, Wilcoxon p = %.4f\n", median(ch$seconds_on_A), median(ct$seconds_on_A), w$p.value))
  cat("\nPseudoreplication: one animal tested 30 times (repeated trials treated as independent) with an individual bias of 0.75 toward A, but no real preference in the population\n")
  fp <- mean(replicate(2000, { bias <- runif(1, 0.3, 0.9); binom.test(rbinom(1, 30, bias), 30, 0.5)$p.value < 0.05 })); cat(sprintf("false positive rate when 30 trials come from ONE random animal: %.2f (should be 0.05)\n", fp))
  fp2 <- mean(replicate(2000, binom.test(rbinom(1, 30, 0.5), 30, 0.5)$p.value < 0.05)); cat(sprintf("false positive rate with 30 independent animals and no real preference: %.2f\n", fp2))
  cat("\nPower: chance of detecting a true 70% preference at p < 0.05 by number of animals (2000 simulations)\n")
  for (n in c(10, 20, 30, 50, 80)) cat(sprintf("n = %2d  power = %.2f\n", n, mean(replicate(2000, binom.test(rbinom(1, n, 0.7), n, 0.5)$p.value < 0.05))))
}
