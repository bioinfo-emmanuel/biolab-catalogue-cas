# ACH-05  Beer's law with an improvised smartphone photometer. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript ach05_photometer.R [standards.csv unknown.csv]
# standards.csv: rel_conc, channel  (mean pixel value, 0 to 255, of the photo of each standard; include the blank as rel_conc 0)
# unknown.csv: channel  (one or more photos of the unknown). Use the colour channel the dye absorbs (red dye: green channel).
set.seed(505)
absorbance <- function(I, I0) -log10(I / I0)

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 2) { st <- read.csv(args[1]); un <- read.csv(args[2]) } else {
    eps <- 0.9; I0 <- 205; cam <- function(A) I0 * 10^(-A) + rnorm(length(A), 0, 0.8)               # pixel value; the camera clips near black
    conc <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5); st <- data.frame(rel_conc = conc, channel = round(pmax(cam(eps * conc + 0.02 * conc^2 * 6) , 4), 1))
    un <- data.frame(channel = round(cam(rep(eps * 0.5 + 0.02 * 0.25 * 6, 3)), 1)); cat("Practice data: 8 standards (relative units 0 to 1.5), an unknown made at 0.50, a slightly curved response\n\n") }
  I0 <- st$channel[st$rel_conc == 0][1]; st$A <- absorbance(st$channel, I0); print(round(st, 3), row.names = FALSE)
  full <- lm(A ~ rel_conc, st); lin <- lm(A ~ rel_conc, st[st$A <= 0.7, ])
  cat(sprintf("\nAll standards: slope %.3f, R squared %.4f.  Standards with A <= 0.7: slope %.3f, R squared %.4f (n = %d)\n", coef(full)[2], summary(full)$r.squared, coef(lin)[2], summary(lin)$r.squared, sum(st$A <= 0.7)))
  Au <- mean(absorbance(un$channel, I0)); x0 <- (Au - coef(lin)[1]) / coef(lin)[2]
  se <- sigma(lin) * sqrt(1 / nrow(un) + 1 / nrow(lin$model) + (x0 - mean(lin$model$rel_conc))^2 / sum((lin$model$rel_conc - mean(lin$model$rel_conc))^2)) / coef(lin)[2]
  cat(sprintf("Unknown: mean absorbance %.3f, relative concentration %.3f, approximate 95%% interval %.3f to %.3f\n", Au, x0, x0 - qt(.975, nrow(lin$model) - 2) * se, x0 + qt(.975, nrow(lin$model) - 2) * se))
  cat("\nLargest absorbance measurable when the darkest pixel value you can trust is 10 (blank", I0, "):", round(absorbance(10, I0), 2), "\n")
  cat("Blank-to-blank pixel scatter is small, but stray light and camera exposure changes are the real limits; repeat each photo three times and compare.\n")
}
