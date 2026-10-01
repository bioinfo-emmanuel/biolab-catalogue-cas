# BOT-06  Testing a dichotomous key: identification accuracy and scorer agreement. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript bot06_key_test.R [tests.csv]
# tests.csv: key, tester, true_species, identified_as   (one row per unknown specimen identified with a key by a tester who did not write it)
set.seed(6)
cohen_kappa <- function(a, b) { lv <- union(a, b); a <- factor(a, lv); b <- factor(b, lv); tab <- table(a, b); n <- sum(tab); po <- sum(diag(tab)) / n; pe <- sum(rowSums(tab) * colSums(tab)) / n^2; (po - pe) / (1 - pe) }
acc_ci <- function(x, n) { b <- binom.test(x, n)$conf.int; c(prop = x / n, lower = b[1], upper = b[2]) }

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) d <- read.csv(args[1], stringsAsFactors = FALSE) else {
    sp <- paste0("species_", LETTERS[1:8])
    mk <- function(key, p_ok, n, confuse) { true <- sample(sp, n, replace = TRUE); ok <- runif(n) < p_ok; id <- true
      for (i in which(!ok)) id[i] <- if (true[i] %in% names(confuse) && runif(1) < 0.7) confuse[[true[i]]] else sample(setdiff(sp, true[i]), 1)
      data.frame(key = key, tester = paste0("tester", (seq_len(n) %% 4) + 1), true_species = true, identified_as = id, stringsAsFactors = FALSE) }
    d <- rbind(mk("key_1", 0.92, 40, list()), mk("key_2", 0.72, 40, list(species_C = "species_D", species_D = "species_C")), mk("key_3", 0.50, 40, list(species_F = "species_G", species_G = "species_F", species_A = "species_B")))
    cat("Practice data: three invented keys for 8 species, each tested with 40 identifications by testers who did not write the key\n\n") }
  d$ok <- d$true_species == d$identified_as
  res <- do.call(rbind, lapply(split(d, d$key), function(x) { r <- acc_ci(sum(x$ok), nrow(x)); data.frame(key = x$key[1], n = nrow(x), correct = sum(x$ok), accuracy = r[1], ci_low = r[2], ci_high = r[3]) })); print(format(res, digits = 3), row.names = FALSE)
  cat("\nComparison of keys (chi-square test of accuracy):\n"); tb <- table(d$key, d$ok); print(tb); print(chisq.test(tb, correct = FALSE))
  cat("\nWhich species fail most often (all keys):\n"); f <- aggregate(ok ~ true_species, d, function(x) c(n = length(x), errors = sum(!x))); f <- do.call(data.frame, f); names(f)[2:3] <- c("n", "errors"); print(f[order(-f$errors), ][1:4, ], row.names = FALSE)
  cat("\nMost common wrong identifications (true -> chosen):\n"); w <- d[!d$ok, ]; print(head(sort(table(paste(w$true_species, "->", w$identified_as)), decreasing = TRUE), 4))
  cat("\nScorer agreement on a character (two people score the same 20 specimens as 'hairy' or 'smooth'): kappa =",
      round(cohen_kappa(c(rep("hairy", 9), rep("smooth", 11)), c(rep("hairy", 7), "smooth", "smooth", rep("smooth", 10), "hairy")), 2), "\n")
  cat("A 95% interval from 40 identifications is still wide: with 30 of 40 correct the interval is", paste(round(acc_ci(30, 40)[2:3], 2), collapse = " to "), "\n")
}
