# BOT-03  Photosynthetic pigments by paper chromatography. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript bot03_pigment_rf.R [bands.csv]
# bands.csv: leaf, run, band (colour or name), spot_cm, front_cm, darkness (optional pixel darkness 0-255 of the band, higher = darker)
# Expected order from the top of the sheet (highest Rf) to the bottom: carotene (yellow-orange), xanthophylls (yellow), chlorophyll a (blue-green), chlorophyll b (yellow-green).
set.seed(3)
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) d <- read.csv(args[1]) else {
    rf0 <- c(carotene = 0.95, xanthophyll = 0.60, chlorophyll_a = 0.40, chlorophyll_b = 0.30)
    mk <- function(leaf, run, ratio) { fr <- runif(1, 8.5, 9.5); do.call(rbind, lapply(names(rf0), function(b) data.frame(leaf = leaf, run = run, band = b, spot_cm = round(rf0[[b]] * fr + rnorm(1, 0, 0.12), 2), front_cm = round(fr, 2),
      darkness = round(c(carotene = 60, xanthophyll = 80, chlorophyll_a = 170, chlorophyll_b = 170 / ratio)[[b]] + rnorm(1, 0, 6), 0)))) }
    d <- rbind(do.call(rbind, lapply(1:4, function(r) mk("sun_leaf", r, 3.2))), do.call(rbind, lapply(1:4, function(r) mk("shade_leaf", r, 2.4))))
    cat("Practice data: sun and shade leaves, 4 runs each, four bands. Rf values and darkness are invented; measure your own on your own paper and solvent.\n\n") }
  d$rf <- d$spot_cm / d$front_cm
  tab <- aggregate(rf ~ leaf + band, d, function(x) c(mean = mean(x), sd = sd(x))); tab <- do.call(data.frame, tab); names(tab)[3:4] <- c("rf_mean", "rf_sd")
  tab <- tab[order(tab$leaf, -tab$rf_mean), ]; print(format(tab, digits = 3), row.names = FALSE)
  if ("darkness" %in% names(d)) {
    w <- reshape(d[d$band %in% c("chlorophyll_a", "chlorophyll_b"), c("leaf", "run", "band", "darkness")], idvar = c("leaf", "run"), timevar = "band", direction = "wide")
    w$ratio_ab <- w$darkness.chlorophyll_a / w$darkness.chlorophyll_b
    cat("\nRatio of band darkness, chlorophyll a to b (a rough index, not a measured concentration ratio):\n"); ag <- do.call(data.frame, aggregate(ratio_ab ~ leaf, w, function(x) c(mean = mean(x), sd = sd(x)))); ag[, -1] <- round(ag[, -1], 2); print(ag, row.names = FALSE)
    tt <- t.test(ratio_ab ~ leaf, w); cat(sprintf("Welch t test between leaf types: p = %.4f\n", tt$p.value)) }
  cat("\nRf order check (highest to lowest) for the first leaf type:\n"); x <- tab[tab$leaf == tab$leaf[1], ]; print(paste(x$band[order(-x$rf_mean)], collapse = " > "))
  cat("Darkness on a photograph is not linear in pigment amount and depends on the camera; use it to compare bands within a sheet.\n")
}
