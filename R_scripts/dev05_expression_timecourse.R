# DEV-05  Gene expression across developmental stages: clustering temporal profiles. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript dev05_expression_timecourse.R [expr.csv geneset.txt]
# expr.csv: first column gene id, then one column per stage in developmental order (normalized expression, not logs).
# geneset.txt: one gene id per line (the category you want to test for enrichment).
# With no arguments the script simulates 900 genes with early, transient, late and flat profiles across 8 stages.
set.seed(325)

simulate_timecourse <- function(n_per = 200, stages = 8, n_flat = 300) {
  t <- seq_len(stages); prof <- list(early = 6 * exp(-(t - 1) / 1.5), transient = 6 * exp(-((t - 4.5)^2) / 2), late = 6 * exp(-(stages - t) / 1.5))
  m <- do.call(rbind, lapply(names(prof), function(nm) t(replicate(n_per, 2^(prof[[nm]] + rnorm(stages, 0, 0.6) + rnorm(1, 2, 0.5))))))
  flat <- t(replicate(n_flat, 2^(rnorm(1, 5, 1.5) + rnorm(stages, 0, 0.5)))); m <- rbind(m, flat)
  rownames(m) <- paste0("g", seq_len(nrow(m))); colnames(m) <- paste0("stage", t)
  list(expr = m, truth = c(rep(names(prof), each = n_per), rep("flat", n_flat)))
}
cluster_profiles <- function(expr, k = 4, min_max = 8) {
  keep <- apply(expr, 1, max) >= min_max; lx <- log2(expr[keep, ] + 1); z <- t(scale(t(lx)))   # standardize each gene
  z <- z[complete.cases(z), ]; hc <- hclust(dist(z), "ward.D2"); list(z = z, cl = cutree(hc, k)) }
enrichment <- function(cl, geneset) { genes <- names(cl); out <- lapply(sort(unique(cl)), function(k) {
    a <- sum(cl == k & genes %in% geneset); b <- sum(cl == k & !genes %in% geneset); c <- sum(cl != k & genes %in% geneset); d <- sum(cl != k & !genes %in% geneset)
    f <- fisher.test(matrix(c(a, b, c, d), 2), alternative = "greater"); data.frame(cluster = k, in_set = a, cluster_size = a + b, odds = unname(f$estimate), p = f$p.value) })
  r <- do.call(rbind, out); r$p_adj <- p.adjust(r$p, "BH"); r }

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE); truth <- NULL
  if (length(args) >= 2) { expr <- as.matrix(read.csv(args[1], row.names = 1)); geneset <- readLines(args[2]) } else {
    s <- simulate_timecourse(); expr <- s$expr; truth <- setNames(s$truth, rownames(expr))
    geneset <- names(truth)[truth == "transient"][1:120]; geneset <- c(geneset, sample(names(truth)[truth != "transient"], 30))   # a category concentrated in transient genes
    cat("Practice data: 900 simulated genes over 8 stages; the gene set has 150 genes, 120 of them transient\n") }
  r <- cluster_profiles(expr, 4); cat(sprintf("\nGenes clustered after filtering: %d\n", nrow(r$z)))
  cat("Cluster sizes:", table(r$cl), "\n\nMean standardized profile per cluster (stages in columns)\n"); print(round(t(sapply(sort(unique(r$cl)), function(k) colMeans(r$z[r$cl == k, , drop = FALSE]))), 2))
  if (!is.null(truth)) { cat("\nAgreement with the simulated profile types\n"); print(table(cluster = r$cl, truth = truth[names(r$cl)])) }
  cat("\nEnrichment of the gene set in each cluster (one-sided Fisher test, BH adjusted)\n"); print(enrichment(r$cl, geneset), digits = 3, row.names = FALSE)
  cat("\nStability: two clusterings on independent random 80% subsets of genes\n")
  a <- sample(rownames(r$z), 0.8 * nrow(r$z)); b <- sample(rownames(r$z), 0.8 * nrow(r$z)); ca <- cluster_profiles(expr[a, ], 4)$cl; cb <- cluster_profiles(expr[b, ], 4)$cl
  sh <- intersect(names(ca), names(cb)); tab <- table(ca[sh], cb[sh]); cat(sprintf("share of shared genes that stay together (best-matching cluster pairs): %.2f\n", sum(apply(tab, 1, max)) / sum(tab)))
  png("profiles.png", 800, 600); par(mfrow = c(2, 2)); for (k in 1:4) matplot(t(r$z[r$cl == k, ]), type = "l", lty = 1, col = "#00000030", main = paste("Cluster", k), xlab = "Stage", ylab = "z-score"); dev.off()
}
