# ENV-06  Terrain analysis from an elevation grid: profile, slope, and exposure. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript env06_terrain.R [grid.csv cellsize_m]
# grid.csv is a matrix of elevations in metres with no header (rows x columns). With no argument the script uses
# R's built-in `volcano` matrix (87 x 61 cells, 10 m spacing; Maunga Whau volcano, Auckland).

terrain_slope <- function(z, cell) {           # central differences; edge cells are dropped
  nr <- nrow(z); nc <- ncol(z)
  dzdx <- (z[2:(nr - 1), 3:nc] - z[2:(nr - 1), 1:(nc - 2)]) / (2 * cell)
  dzdy <- (z[3:nr, 2:(nc - 1)] - z[1:(nr - 2), 2:(nc - 1)]) / (2 * cell)
  atan(sqrt(dzdx^2 + dzdy^2)) * 180 / pi
}
profile_line <- function(z, cell, row) data.frame(distance_m = (seq_len(ncol(z)) - 1) * cell, elevation_m = z[row, ])
aggregate_grid <- function(z, k) { nr <- (nrow(z) %/% k) * k; nc <- (ncol(z) %/% k) * k
  m <- z[1:nr, 1:nc]; out <- matrix(NA, nr / k, nc / k)
  for (i in 1:(nr / k)) for (j in 1:(nc / k)) out[i, j] <- mean(m[((i - 1) * k + 1):(i * k), ((j - 1) * k + 1):(j * k)]); out }

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 2) { z <- as.matrix(read.csv(args[1], header = FALSE)); cell <- as.numeric(args[2]) } else { z <- datasets::volcano; cell <- 10 }
  cat(sprintf("Grid: %d rows x %d columns, cell size %g m, elevation %g to %g m\n", nrow(z), ncol(z), cell, min(z), max(z)))
  s <- terrain_slope(z, cell)
  cat(sprintf("\nSlope (degrees): mean %.1f, median %.1f, maximum %.1f\n", mean(s), median(s), max(s)))
  for (th in c(15, 25, 35)) cat(sprintf("Area with slope above %d degrees: %.1f%% of cells (%.1f hectares)\n", th, 100 * mean(s > th), sum(s > th) * cell^2 / 1e4))
  zi <- z[2:(nrow(z) - 1), 2:(ncol(z) - 1)]
  for (lim in c(110, 130)) cat(sprintf("Cells below %d m elevation: %.1f%%; both below %d m and slope under 5 degrees (flood-prone flat land): %.1f%%\n", lim, 100 * mean(zi < lim), lim, 100 * mean(zi < lim & s < 5)))
  p <- profile_line(z, cell, 44); cat(sprintf("\nTransect along row 44: %d points over %d m, rise %.0f m to the highest point, highest point %.0f m at %d m\n",
      nrow(p), max(p$distance_m), max(p$elevation_m) - min(p$elevation_m), max(p$elevation_m), p$distance_m[which.max(p$elevation_m)]))
  cat("\nScale dependence: the same terrain at coarser resolution\n")
  for (k in c(1, 2, 4)) { zz <- if (k == 1) z else aggregate_grid(z, k); ss <- terrain_slope(zz, cell * k)
    cat(sprintf("cell size %2d m: mean slope %.1f degrees, area above 25 degrees %.1f%%\n", cell * k, mean(ss), 100 * mean(ss > 25))) }
  png("terrain_contours.png", 700, 500); contour(z, main = "Elevation (m)"); dev.off()
}
