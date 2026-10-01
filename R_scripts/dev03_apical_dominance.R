# DEV-03  Apical dominance and axillary bud outgrowth. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript dev03_apical_dominance.R [buds.csv]
# buds.csv: plant, treatment (intact, decapitated, decapitated_auxin), day, longest_lateral_mm, laterals_over_5mm
# Assign plants to treatments at random (see random_assign) before cutting anything.
set.seed(3)
random_assign <- function(plants, treatments) { n <- length(plants); data.frame(plant = plants, treatment = sample(rep(treatments, length.out = n))) }
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) d <- read.csv(args[1]) else {
    a <- random_assign(paste0("P", 1:30), c("intact", "decapitated", "decapitated_auxin")); days <- c(0, 2, 4, 6, 8, 10); g <- c(intact = 0.12, decapitated = 1.6, decapitated_auxin = 0.35)
    d <- do.call(rbind, lapply(seq_len(nrow(a)), function(i) { r <- g[[a$treatment[i]]] * rnorm(1, 1, 0.25); L <- pmax(0, r * days * (1 + 0.05 * days) * rnorm(length(days), 1, 0.1)); data.frame(plant = a$plant[i], treatment = a$treatment[i], day = days, longest_lateral_mm = round(L, 1), laterals_over_5mm = pmin(4, floor(L / 5))) }))
    cat("Practice data: 30 invented plants randomly assigned to three treatments (10 each), measured every 2 days for 10 days. Growth rates are invented.\n\n") }
  last <- d[d$day == max(d$day), ]; cat("Longest lateral shoot on the last day (mm), by treatment:\n")
  res <- do.call(rbind, lapply(split(last, last$treatment), function(x) data.frame(treatment = x$treatment[1], n = nrow(x), mean = mean(x$longest_lateral_mm), sd = sd(x$longest_lateral_mm), median = median(x$longest_lateral_mm), share_with_lateral_over_5mm = mean(x$laterals_over_5mm > 0)))); print(format(res, digits = 3), row.names = FALSE)
  k <- kruskal.test(longest_lateral_mm ~ treatment, last); cat(sprintf("\nKruskal-Wallis test across treatments: chi-square = %.1f, p = %.2e\n", k$statistic, k$p.value))
  cat("Pairwise Wilcoxon tests (Holm-adjusted p):\n"); print(signif(pairwise.wilcox.test(last$longest_lateral_mm, last$treatment, p.adjust.method = "holm", exact = FALSE)$p.value, 2))
  fp <- fisher.test(table(last$treatment[last$treatment != "decapitated_auxin"], last$laterals_over_5mm[last$treatment != "decapitated_auxin"] > 0)); cat(sprintf("Intact vs decapitated, share of plants with a lateral over 5 mm: Fisher exact p = %.4f\n", fp$p.value))
  cat("\nGrowth rate of the longest lateral (mm per day, slope over days), mean by treatment:\n"); sl <- sapply(split(d, d$plant), function(x) coef(lm(longest_lateral_mm ~ day, x))[2]); tr <- tapply(d$treatment, d$plant, `[`, 1); print(round(tapply(sl, tr, mean), 2))
  cat("\nRandomization check: the first six plants were assigned to ", paste(random_assign(paste0("P", 1:6), c("intact", "decapitated", "decapitated_auxin"))$treatment, collapse = ", "), "\n")
}
