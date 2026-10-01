# ZOO-01  Comparative invertebrate body plans: character matrix, distances, and scorer agreement. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript zoo01_body_plan_matrix.R [matrix.csv]
# matrix.csv: first column taxon, then one column per character scored 0 (absent) or 1 (present).
# The built-in matrix scores textbook states for five common animals and is a starting point for checking student scoring, not a substitute for it.
set.seed(1)
demo <- data.frame(taxon = c("earthworm", "shrimp", "snail", "mussel", "squid"),
  bilateral_symmetry = c(1, 1, 1, 1, 1), external_segmentation = c(1, 1, 0, 0, 0), jointed_appendages = c(0, 1, 0, 0, 0), chitinous_exoskeleton = c(0, 1, 0, 0, 0),
  calcareous_shell = c(0, 0, 1, 1, 0), muscular_foot = c(0, 0, 1, 1, 0), mantle = c(0, 0, 1, 1, 1), radula = c(0, 0, 1, 0, 1),
  closed_circulatory_system = c(1, 0, 0, 0, 1), extensive_true_coelom = c(1, 0, 0, 0, 0), ventral_nerve_cord = c(1, 1, 0, 0, 0), complete_digestive_tract = c(1, 1, 1, 1, 1))
simple_match <- function(a, b) mean(a == b)
cohen_kappa <- function(a, b) { tab <- table(factor(a, 0:1), factor(b, 0:1)); n <- sum(tab); po <- sum(diag(tab)) / n; pe <- sum(rowSums(tab) * colSums(tab)) / n^2; if (pe == 1) 1 else (po - pe) / (1 - pe) }
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE); m <- if (length(args) >= 1) read.csv(args[1]) else demo
  X <- as.matrix(m[, -1]); rownames(X) <- m$taxon
  cat(sprintf("%d taxa x %d characters\n\nCharacter matrix:\n", nrow(X), ncol(X))); print(X)
  D <- as.matrix(dist(X, method = "manhattan")); cat("\nNumber of characters that differ between each pair of taxa:\n"); print(D)
  cat("\nCharacters shared as present by exactly two taxa (candidates for a shared history or convergence, to be argued by the students):\n")
  for (j in seq_len(ncol(X))) if (sum(X[, j]) == 2) cat(sprintf("%-28s %s\n", colnames(X)[j], paste(rownames(X)[X[, j] == 1], collapse = " + ")))
  cat("\nCharacters that are present in every taxon (uninformative for grouping):", paste(colnames(X)[colSums(X) == nrow(X)], collapse = ", "), "\n")
  cat("Characters unique to one taxon:", paste(colnames(X)[colSums(X) == 1], collapse = ", "), "\n")
  hc <- hclust(dist(X, method = "manhattan"), "average"); cat("\nAverage-linkage clustering order:", paste(rownames(X)[hc$order], collapse = " - "), "\n")
  cat("Merges (height = characters differing):\n"); print(data.frame(merge_step = seq_len(nrow(hc$merge)), height = hc$height))
  cat("\nWarning: a distance tree from 12 characters and 5 taxa shows similarity, which mixes shared ancestry with convergence.\n")
  cat("\nObserver agreement: two students score the same 12 characters for one animal (10 agree):\n")
  a <- demo[1, -1] |> unlist() |> as.integer(); b <- a; b[c(4, 9)] <- 1 - b[c(4, 9)]; cat(sprintf("simple agreement %.2f, kappa %.2f\n", simple_match(a, b), cohen_kappa(a, b)))
}
