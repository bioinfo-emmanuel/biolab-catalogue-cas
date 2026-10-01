# AP-03  Mammalian kidney: dissection measurements and comparison with histology. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript ap03_kidney_measurements.R [kidneys.csv]
# kidneys.csv: kidney, species, observer, mass_g, length_mm, cortex_mm, medulla_mm, pyramids, glomerulus_um (optional; from the histology slide, calibrated)
# Cortex and medulla thickness are measured on a cut face at the same place (mid-height of a pyramid), with a ruler or calliper.
set.seed(3)
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) d <- read.csv(args[1]) else {
    sp <- rep(c("goat", "pig"), each = 6); mk <- function(obs, sh) data.frame(kidney = paste0("K", 1:12), species = sp, observer = obs, mass_g = round(ifelse(sp == "goat", rnorm(12, 75, 8), rnorm(12, 140, 12)), 0),
      length_mm = round(ifelse(sp == "goat", rnorm(12, 68, 4), rnorm(12, 105, 6)), 0), cortex_mm = round(ifelse(sp == "goat", rnorm(12, 8, 0.8), rnorm(12, 11, 1)) + sh, 1), medulla_mm = round(ifelse(sp == "goat", rnorm(12, 16, 1.5), rnorm(12, 22, 2)) + sh, 1), pyramids = ifelse(sp == "goat", 1, sample(7:10, 12, TRUE)))
    a <- mk("A", 0); b <- a; b$observer <- "B"; b$cortex_mm <- round(a$cortex_mm + rnorm(12, 0, 0.5), 1); b$medulla_mm <- round(a$medulla_mm + rnorm(12, 0, 0.7), 1); d <- rbind(a, b)
    cat("Practice data: 6 invented goat-like and 6 invented pig-like kidneys measured by two observers. Sizes and pyramid counts are invented for practice.\n\n") }
  d$cm_ratio <- d$cortex_mm / d$medulla_mm; a <- d[d$observer == unique(d$observer)[1], ]
  res <- do.call(rbind, lapply(split(a, a$species), function(x) data.frame(species = x$species[1], n = nrow(x), mass_g = mean(x$mass_g), length_mm = mean(x$length_mm), cortex_mm = mean(x$cortex_mm), medulla_mm = mean(x$medulla_mm), cortex_to_medulla = mean(x$cm_ratio), pyramids = median(x$pyramids)))); print(format(res, digits = 3), row.names = FALSE)
  w <- wilcox.test(cm_ratio ~ species, a, exact = FALSE); cat(sprintf("\nCortex:medulla ratio, species compared: Wilcoxon p = %.4f (n = %d per species)\n", w$p.value, sum(a$species == a$species[1])))
  cat("Mass against length across all kidneys (two species pooled, so this mostly shows species size):"); f <- lm(log(mass_g) ~ log(length_mm), a); cat(sprintf(" log-log slope %.2f, R squared %.2f; geometric similarity would give a slope of 3\n", coef(f)[2], summary(f)$r.squared))
  o <- reshape(d[, c("kidney", "observer", "cortex_mm")], idvar = "kidney", timevar = "observer", direction = "wide"); df <- o[[2]] - o[[3]]
  cat(sprintf("\nObserver agreement on cortex thickness: mean difference %.2f mm, limits of agreement %.2f to %.2f mm; within-kidney sd of the two readings about %.2f mm\n", mean(df), mean(df) - 1.96 * sd(df), mean(df) + 1.96 * sd(df), sd(df) / sqrt(2)))
  cat("A single cut face shows only one plane of a three-dimensional organ; report where the cut was made.\n")
}
