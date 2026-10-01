# BOT-05  Plant anatomy by hand sectioning: monocot versus dicot stems. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript bot05_vascular_anatomy.R [sections.csv vessels.csv]
# sections.csv: plant, group (monocot or dicot), section, bundles, field_diameter_mm (calibrated diameter of the field of view used to count)
# vessels.csv: plant, group, vessel_um  (widest diameter of individual xylem vessels, calibrated)
# Calibrate the field once with a millimetre ruler on the stage at the same objective, and record the field diameter in mm.
set.seed(5)
field_area <- function(d_mm) pi * (d_mm / 2)^2
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 2) { s <- read.csv(args[1]); v <- read.csv(args[2]) } else {
    mk <- function(plant, grp, dens, n) data.frame(plant = plant, group = grp, section = 1:n, field_diameter_mm = 1.8, bundles = rpois(n, dens * field_area(1.8)))
    s <- rbind(mk("plant_M1", "monocot", 9, 10), mk("plant_M2", "monocot", 12, 10), mk("plant_D1", "dicot", 3, 10), mk("plant_D2", "dicot", 4, 10))
    v <- rbind(data.frame(plant = "plant_M1", group = "monocot", vessel_um = round(rlnorm(30, log(85), 0.25))), data.frame(plant = "plant_D1", group = "dicot", vessel_um = round(rlnorm(30, log(48), 0.25))))
    cat("Practice data: bundles per field at 1.8 mm field diameter for two invented monocot and two invented dicot stems, and 30 vessel widths per group\n\n") }
  s$density <- s$bundles / field_area(s$field_diameter_mm)                       # bundles per mm2
  res <- do.call(rbind, lapply(split(s, s$group), function(x) data.frame(group = x$group[1], sections = nrow(x), mean_per_mm2 = mean(x$density), sd = sd(x$density))))
  print(format(res, digits = 3), row.names = FALSE)
  w <- wilcox.test(density ~ group, s, exact = FALSE); cat(sprintf("\nBundle density, monocot vs dicot: Wilcoxon p = %.5f\n", w$p.value))
  cat("\nPlant means (bundles per mm2), the independent units:\n"); pm <- aggregate(density ~ plant + group, s, mean); pm$density <- round(pm$density, 2); print(pm, row.names = FALSE)
  cat("Note: sections from the same stem are not independent; the statistically sound unit is the plant, so the test above overstates the evidence when only two plants per group are used.\n")
  vr <- do.call(rbind, lapply(split(v, v$group), function(x) data.frame(group = x$group[1], n = nrow(x), median_um = median(x$vessel_um), iqr = IQR(x$vessel_um))))
  cat("\nVessel diameters (um):\n"); print(vr, row.names = FALSE); wv <- wilcox.test(vessel_um ~ group, v, exact = FALSE); cat(sprintf("Wilcoxon p = %.2e\n", wv$p.value))
  cat("\nField diameter check: a 2.0 mm field has area", round(field_area(2.0), 2), "mm2, so 10 bundles in that field is", round(10 / field_area(2.0), 1), "per mm2\n")
}
