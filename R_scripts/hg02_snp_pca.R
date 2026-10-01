# HG-02  Population structure: PCA of SNP genotypes and assignment of unknowns. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript hg02_snp_pca.R [genotypes.csv]
# genotypes.csv: rows are individuals, columns are SNPs coded 0, 1, 2 (copies of one allele), plus a first column
# named population (use NA for unknown individuals). With no argument the script simulates three populations
# (Balding-Nichols model) so it runs offline. Never use students' own DNA.
set.seed(2025)

simulate_snps <- function(n_snp = 500, n_per = 40, fst = 0.05, n_unknown = 6) {
  anc <- runif(n_snp, 0.1, 0.9); pops <- c("PopA", "PopB", "PopC")
  fr <- sapply(pops, function(p) rbeta(n_snp, anc * (1 - fst) / fst, (1 - anc) * (1 - fst) / fst))
  gen <- function(pop, n) t(sapply(seq_len(n), function(i) rbinom(n_snp, 2, fr[, pop])))
  known <- do.call(rbind, lapply(pops, function(p) gen(p, n_per))); unk_true <- rep(pops[2], n_unknown); unk <- gen(pops[2], n_unknown)
  g <- rbind(known, unk); colnames(g) <- paste0("snp", seq_len(n_snp))
  data.frame(population = c(rep(pops, each = n_per), rep(NA, n_unknown)), true_pop = c(rep(pops, each = n_per), unk_true), g, check.names = FALSE) }

pca_snps <- function(G) { keep <- apply(G, 2, function(x) sd(x) > 0); prcomp(G[, keep], center = TRUE, scale. = TRUE) }
assign_nearest <- function(pc, pop, k = 2) { known <- !is.na(pop); cen <- apply(pc[known, 1:k, drop = FALSE], 2, function(x) tapply(x, pop[known], mean))
  apply(pc[!known, 1:k, drop = FALSE], 1, function(u) rownames(cen)[which.min(colSums((t(cen) - u)^2))]) }

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) { d <- read.csv(args[1], check.names = FALSE); d$true_pop <- NA } else { d <- simulate_snps(); cat("Practice data: 3 simulated populations (Fst about 0.05), 40 individuals each, 500 SNPs, 6 unknowns\n") }
  G <- as.matrix(d[, setdiff(names(d), c("population", "true_pop"))]); pop <- d$population
  p <- pca_snps(G); ve <- 100 * p$sdev^2 / sum(p$sdev^2)
  cat(sprintf("\nVariance explained: PC1 %.1f%%, PC2 %.1f%%, PC3 %.1f%%, PC4 %.1f%%\n", ve[1], ve[2], ve[3], ve[4]))
  cat("\nMean PC1 and PC2 by population:\n"); m <- aggregate(as.data.frame(p$x[!is.na(pop), 1:2]), list(population = pop[!is.na(pop)]), mean); m[, -1] <- round(m[, -1], 2); print(m, row.names = FALSE)
  a <- assign_nearest(p$x, pop); cat("\nAssignments of the unknown individuals:", a, "\n"); if (!all(is.na(d$true_pop[is.na(pop)]))) cat("Truth (simulation only):", d$true_pop[is.na(pop)], "\n")
  cat("\nDoes the number of SNPs matter? Assignment accuracy on 60 held-out individuals (30 repeats each)\n")
  for (ns in c(20, 50, 200, 500)) { acc <- mean(replicate(30, { s <- simulate_snps(n_snp = ns, n_unknown = 60); Gs <- as.matrix(s[, setdiff(names(s), c("population", "true_pop"))]);
      pc <- pca_snps(Gs)$x; aa <- assign_nearest(pc, s$population); mean(aa == s$true_pop[is.na(s$population)]) })); cat(sprintf("%3d SNPs: accuracy %.2f\n", ns, acc)) }
  png("snp_pca.png", 600, 500); plot(p$x[, 1], p$x[, 2], col = as.integer(factor(pop)) + 1, pch = 19, xlab = "PC1", ylab = "PC2"); dev.off()
}
