# ZOO-05  Animal tissue types and cell measurements from prepared slides. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript zoo05_cell_measurements.R [cells.csv]
# cells.csv: tissue, cell, observer, length, width   (long axis and short axis of the cell or nucleus as cut in the section; any consistent unit; use micrometres if the field was calibrated, otherwise eyepiece divisions or millimetres on a photograph)
# Ratios such as length:width need no calibration. Absolute sizes need the field to be calibrated with a ruler at the same objective.
set.seed(5)
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) d <- read.csv(args[1]) else {
    spec <- data.frame(tissue = c("simple_squamous", "simple_cuboidal", "simple_columnar", "smooth_muscle", "skeletal_muscle"), len = c(30, 15, 25, 70, 400), wid = c(4, 14, 9, 6, 45))
    d <- do.call(rbind, lapply(seq_len(nrow(spec)), function(i) { n <- 30; l <- rlnorm(n, log(spec$len[i]), 0.18); w <- rlnorm(n, log(spec$wid[i]), 0.18)
      rbind(data.frame(tissue = spec$tissue[i], cell = 1:n, observer = "A", length = round(l, 1), width = round(w, 1)), data.frame(tissue = spec$tissue[i], cell = 1:n, observer = "B", length = round(l * rnorm(n, 1, 0.05), 1), width = round(w * rnorm(n, 1, 0.05), 1))) }))
    cat("Practice data: 30 invented cells per tissue for five tissues, each measured by two observers (units invented)\n\n") }
  d$ratio <- d$length / d$width
  a <- d[d$observer == unique(d$observer)[1], ]
  res <- do.call(rbind, lapply(split(a, a$tissue), function(x) data.frame(tissue = x$tissue[1], n = nrow(x), median_length = median(x$length), median_width = median(x$width), median_ratio = median(x$ratio), iqr_ratio = IQR(x$ratio)))); print(format(res, digits = 3), row.names = FALSE)
  k <- kruskal.test(ratio ~ tissue, a); cat(sprintf("\nKruskal-Wallis test of length:width ratio across tissues: chi-square = %.1f, df = %d, p = %.2e\n", k$statistic, k$parameter, k$p.value))
  cat("\nPairwise Wilcoxon tests with Holm adjustment (p values):\n"); print(signif(pairwise.wilcox.test(a$ratio, a$tissue, p.adjust.method = "holm", exact = FALSE)$p.value, 2))
  cat("\nObserver agreement on the same cells (ratio): Bland-Altman style\n"); w <- reshape(d[, c("tissue", "cell", "observer", "ratio")], idvar = c("tissue", "cell"), timevar = "observer", direction = "wide"); ob <- names(w)[grepl("ratio", names(w))]
  df <- w[[ob[1]]] - w[[ob[2]]]; cat(sprintf("mean difference %.3f, limits of agreement %.2f to %.2f, correlation %.3f\n", mean(df), mean(df) - 1.96 * sd(df), mean(df) + 1.96 * sd(df), cor(w[[ob[1]]], w[[ob[2]]])))
  cat("\nHow many cells are enough? Standard error of the median ratio for the first tissue by sample size (bootstrap):\n"); x0 <- a$ratio[a$tissue == a$tissue[1]]; for (n in c(5, 10, 20, 30)) cat(sprintf("n = %2d  se of median = %.3f\n", n, sd(replicate(500, median(sample(x0, n, replace = TRUE))))))
}
