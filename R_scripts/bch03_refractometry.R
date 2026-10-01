# BCH-03  Sugar content of beverages by refractometry: calibration and analysis. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript bch03_refractometry.R [standards.csv beverages.csv]
# standards.csv: sucrose_pct_w_w, brix_reading, temp_c.   beverages.csv: beverage, brix_reading, temp_c, label_g_per_100ml (optional)
# A Brix scale is defined as percent sucrose by mass, so an ideal refractometer reads the standard's percentage directly.
set.seed(33)
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 2) { st <- read.csv(args[1]); bv <- read.csv(args[2]) } else {
    pct <- c(0, 2, 4, 6, 8, 10, 12, 15); st <- data.frame(sucrose_pct_w_w = pct, brix_reading = round(0.3 + 0.985 * pct + rnorm(length(pct), 0, 0.12), 1), temp_c = 27)
    bv <- data.frame(beverage = rep(c("cola", "juice", "sports_drink", "tea_drink"), each = 3), brix_reading = round(rep(c(10.6, 11.8, 6.2, 8.4), each = 3) + rnorm(12, 0, 0.12), 1), temp_c = 27, label_g_per_100ml = rep(c(10.6, 11.0, 6.0, 8.0), each = 3))
    cat("Practice data: 8 sucrose standards (0 to 15 percent w/w) and four beverages read in triplicate, all at 27 C. Reading offset and slope are simulated.\n\n") }
  fit <- lm(brix_reading ~ sucrose_pct_w_w, st); b <- coef(fit)
  cat(sprintf("Calibration: reading = %.3f + %.3f x sucrose (%% w/w); R squared %.4f; residual sd %.3f\n", b[1], b[2], summary(fit)$r.squared, sigma(fit)))
  cat(sprintf("Zero check: distilled water reads %.2f (should be 0 on a correct instrument)\n", st$brix_reading[st$sucrose_pct_w_w == 0][1]))
  bv$sucrose_pct <- (bv$brix_reading - b[1]) / b[2]
  # convert % w/w to g per 100 mL using the density of a sucrose solution: rho approx 0.9982 + 0.0038 * pct (g/mL), an approximation that is within about 0.6% of tabulated values below 20% w/w
  bv$g_per_100ml <- bv$sucrose_pct * (0.9982 + 0.0038 * bv$sucrose_pct)
  res <- do.call(rbind, lapply(split(bv, bv$beverage), function(x) data.frame(beverage = x$beverage[1], n = nrow(x), pct_w_w = mean(x$sucrose_pct), g_per_100ml = mean(x$g_per_100ml), sd = sd(x$g_per_100ml), label = x$label_g_per_100ml[1])))
  res$diff_from_label <- res$g_per_100ml - res$label; print(format(res, digits = 3), row.names = FALSE)
  cat("\nWhat the reading includes: a refractometer responds to all dissolved solids (sugars, acids, salts), so the result is an upper limit for sugar in a mixed beverage.\n")
  cat("Temperature effect: refractive index changes with temperature. Read standards and samples at the same temperature, or use an instrument with automatic temperature compensation and check its stated range.\n")
  cat("Queue planning: 3 groups x (8 standards + 4 beverages x 3 reads) = ", 3 * (8 + 12), "readings, about 2 minutes each with rinsing, or about", 3 * 20 * 2 / 60, "hours of instrument time in total for one refractometer.\n")
}
