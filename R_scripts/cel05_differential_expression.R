# CEL-05  Differential gene expression from count data. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript cel05_differential_expression.R [counts.csv]
# counts.csv: first column gene id, then one column per sample named <group>_<replicate> (for example ctrl_1, treat_1).
# With no argument the script simulates 3000 genes, 5 samples per group, 10% truly changed, so it runs offline.
# This is a teaching pipeline (log-CPM and a t-test). For real projects use edgeR, DESeq2 or limma.
set.seed(2026)

simulate_counts <- function(n_genes = 3000, n_per = 5, frac_de = 0.1, lfc = 1.5) {
  base <- rlnorm(n_genes, log(10), 1.6); de <- rep(FALSE, n_genes); de[sample(n_genes, n_genes * frac_de)] <- TRUE
  eff <- ifelse(de, sample(c(-1, 1), n_genes, TRUE) * lfc, 0)
  lib <- runif(2 * n_per, 0.7, 1.3); grp <- rep(c("ctrl", "treat"), each = n_per)
  mu <- outer(base, lib) * 2^outer(eff, as.numeric(grp == "treat"))
  cnt <- matrix(rnbinom(length(mu), size = 30, mu = mu), n_genes)
  colnames(cnt) <- paste0(grp, "_", rep(1:n_per, 2)); rownames(cnt) <- paste0("g", seq_len(n_genes)); list(counts = cnt, de = de)
}
run_de <- function(cnt, group) {
  keep <- rowSums(sweep(cnt, 2, colSums(cnt) / 1e6, "/") > 1) >= 3; cnt <- cnt[keep, ]
  lcpm <- log2(sweep(cnt, 2, colSums(cnt) / 1e6, "/") + 1); g1 <- unique(group)[1]; g2 <- unique(group)[2]
  res <- t(apply(lcpm, 1, function(y) { tt <- t.test(y[group == g2], y[group == g1]); c(log2fc = mean(y[group == g2]) - mean(y[group == g1]), p = tt$p.value) }))
  res <- data.frame(gene = rownames(res), res, row.names = NULL); res$padj <- p.adjust(res$p, "BH"); res
}

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE); truth <- NULL
  if (length(args) >= 1) { d <- read.csv(args[1], row.names = 1); cnt <- as.matrix(d) } else { s <- simulate_counts(); cnt <- s$counts; truth <- s$de; names(truth) <- rownames(cnt) }
  group <- sub("_.*", "", colnames(cnt)); cat(sprintf("%d genes, groups: %s\n", nrow(cnt), paste(names(table(group)), table(group), collapse = ", ")))
  res <- run_de(cnt, group); cat(sprintf("Genes kept after the low-count filter: %d\n", nrow(res)))
  cat(sprintf("Raw p < 0.05: %d genes.  BH-adjusted p < 0.05: %d genes.\n", sum(res$p < 0.05), sum(res$padj < 0.05)))
  if (!is.null(truth)) { called <- res$padj < 0.05; t <- truth[res$gene]
    cat(sprintf("Known truth (simulation only): true positives %d, false positives %d, realized false discovery proportion %.2f, sensitivity %.2f\n", sum(called & t), sum(called & !t), sum(called & !t) / max(sum(called), 1), sum(called & t) / sum(t)))
    raw <- res$p < 0.05; cat(sprintf("Without correction the same threshold gives %d false positives among %d calls.\n", sum(raw & !t), sum(raw))) }
  cat("\nLabel shuffling (no true difference): significant genes in 20 random relabellings\n")
  sh <- sapply(1:20, function(i) { g <- sample(group); r <- run_de(cnt, g); c(raw = sum(r$p < 0.05), bh = sum(r$padj < 0.05)) })
  cat(sprintf("raw p < 0.05: mean %.0f genes (about %.0f expected by chance); BH < 0.05: mean %.1f genes\n", mean(sh["raw", ]), 0.05 * nrow(res), mean(sh["bh", ])))
  cat("\nTop genes:\n"); print(head(res[order(res$p), ], 6), digits = 3, row.names = FALSE)
  png("volcano.png", 700, 550); plot(res$log2fc, -log10(res$p), pch = 20, col = ifelse(res$padj < 0.05, "red", "grey60"), xlab = "log2 fold change", ylab = "-log10 p"); dev.off()
}
