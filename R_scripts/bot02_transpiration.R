# BOT-02  Transpiration rate under light, wind, and humidity treatments. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript bot02_transpiration.R [shoots.csv]
# shoots.csv: treatment, shoot, minutes, loss_g, leaf_area_cm2   (mass lost by a leafy shoot in a stoppered flask of water, or water taken up in mL, over `minutes`)
# Leaf area can be found from a paper tracing: leaf_area = (mass of cut-out / mass of a reference square) x reference area (see leaf_area_from_paper).
set.seed(2)
leaf_area_from_paper <- function(cutout_g, ref_square_g, ref_area_cm2) cutout_g / ref_square_g * ref_area_cm2
rate_per_area <- function(d) d$loss_g / (d$minutes / 60) / d$leaf_area_cm2                # g of water per cm2 of leaf per hour (1 g of water is about 1 mL)

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) d <- read.csv(args[1]) else {
    mu <- c(control = 0.0012, fan = 0.0021, lamp = 0.0016, humid_bag = 0.00045); n <- 6
    d <- do.call(rbind, lapply(names(mu), function(tr) { area <- rnorm(n, 380, 60); data.frame(treatment = tr, shoot = 1:n, minutes = 120, leaf_area_cm2 = round(area, 0), loss_g = round(mu[[tr]] * area * 2 * rnorm(n, 1, 0.12), 1)) }))
    cat("Practice data: 4 treatments x 6 shoots, 120 minutes, mass read to 0.1 g. Rates are invented for practice.\n\n") }
  d$rate <- rate_per_area(d)
  res <- do.call(rbind, lapply(split(d, d$treatment), function(x) data.frame(treatment = x$treatment[1], n = nrow(x), mean_rate = mean(x$rate), sd = sd(x$rate), mean_loss_g = mean(x$loss_g))))
  res$mean_rate <- res$mean_rate * 1000; res$sd <- res$sd * 1000; names(res)[3:4] <- c("mean_mg_cm2_h", "sd_mg_cm2_h"); print(format(res, digits = 3), row.names = FALSE)
  d$treatment <- factor(d$treatment, levels = unique(d$treatment)); a <- aov(rate ~ treatment, d); print(summary(a))
  cat("\nPairwise comparisons with Holm adjustment (Welch t tests on the rate per area):\n"); print(pairwise.t.test(d$rate, d$treatment, p.adjust.method = "holm", pool.sd = FALSE))
  cat(sprintf("\nReading error: with a balance reading to 0.1 g and two readings, the mass loss has an uncertainty of about %.0f%% for a 1.0 g loss and %.0f%% for a 3.0 g loss.\n", 100 * sqrt(2) * 0.1 / sqrt(12) / 1, 100 * sqrt(2) * 0.1 / sqrt(12) / 3))
  cat("Why normalize: shoots differ in leaf area (here about", round(100 * sd(d$leaf_area_cm2) / mean(d$leaf_area_cm2)), "% sd), and raw loss per shoot mixes leaf area with treatment.\n")
  cat("Check of the paper method: a 100 cm2 reference square weighs 0.80 g and a leaf cut-out 0.36 g gives", round(leaf_area_from_paper(0.36, 0.80, 100), 1), "cm2\n")
}
