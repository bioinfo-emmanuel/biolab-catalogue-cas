# ZOO-02  Gut length, diet, and body size in market fish. Base R only: no packages to install.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
#
# YOUR DATA: a CSV with one row per fish and these columns
#   species, fish, standard_length_mm, gut_length_mm
#   standard_length_mm  snout tip to the base of the tail fin
#   gut_length_mm       esophagus to anus, straightened without stretching
# The fish is the replicate. Relative gut length = gut length / standard length.
#
# Run in RStudio: set data_file below and click Source. From a terminal: Rscript gut_length_analysis.R mydata.csv
# With no file the script makes PRACTICE data (simulated). Never report them as results.

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1) args[1] else NA

if (is.na(data_file)) {
  message("No data file given: using simulated PRACTICE data.")
  set.seed(6)
  rel <- c(species_A = 2.6, species_B = 1.4, species_C = 0.8)      # built-in relative gut lengths
  d <- do.call(rbind, lapply(names(rel), function(sp) {
    sl <- round(rnorm(12, 180, 25)); data.frame(species = sp, fish = 1:12, standard_length_mm = sl,
                                                gut_length_mm = round(sl * rel[[sp]] * rnorm(12, 1, 0.10)))
  }))
} else {
  d <- read.csv(data_file, stringsAsFactors = FALSE)
}
d$relative <- d$gut_length_mm / d$standard_length_mm
d$species <- factor(d$species)

cat("\n1. Relative gut length by species (mean, sd, n)\n")
print(do.call(rbind, lapply(split(d, d$species), function(x) data.frame(species = x$species[1], n = nrow(x),
      mean = round(mean(x$relative), 2), sd = round(sd(x$relative), 2)))), row.names = FALSE)

cat("\n2. Species differences: Kruskal-Wallis, then pairwise Wilcoxon with Holm adjustment\n")
kw <- kruskal.test(relative ~ species, data = d)
cat(sprintf("Kruskal-Wallis chi-squared = %.1f, df = %d, p = %.2e\n", kw$statistic, kw$parameter, kw$p.value))
print(signif(pairwise.wilcox.test(d$relative, d$species, p.adjust.method = "holm", exact = FALSE)$p.value, 2))

cat("\n3. Within each species: Spearman correlation of relative gut length with standard length\n")
print(do.call(rbind, lapply(split(d, d$species), function(x) { s <- suppressWarnings(cor.test(x$standard_length_mm, x$relative, method = "spearman", exact = FALSE))
      data.frame(species = x$species[1], rho = round(unname(s$estimate), 2), p = round(s$p.value, 3)) })), row.names = FALSE)

cat("\n4. Extension: log gut length on log standard length with species as a factor\n")
m <- lm(log(gut_length_mm) ~ log(standard_length_mm) + species, data = d)
print(round(summary(m)$coefficients, 3))
cat("A slope near 1 means gut length grows in proportion to body length, so relative gut length is a fair index.\n")

boxplot(relative ~ species, data = d, ylab = "Gut length / standard length", main = "Relative gut length")
