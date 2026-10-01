# ACH-04  Conductivity calibration, detection limit, and standard addition. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript ach04_conductivity_calibration.R [standards.csv blanks.csv unknown.csv addition.csv]
# standards.csv: conc_mg_l, cond_us_cm, temp_c.  blanks.csv: cond_us_cm.  unknown.csv: cond_us_cm, temp_c.
# addition.csv: added_mg_l, cond_us_cm, temp_c  (standard addition on a natural water sample).
set.seed(404)
to25 <- function(k, t, alpha = 0.02) k / (1 + alpha * (t - 25))      # rough temperature compensation, about 2% per degree C

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 4) { st <- read.csv(args[1]); bl <- read.csv(args[2]); un <- read.csv(args[3]); ad <- read.csv(args[4]) } else {
    slope <- 2.0; T0 <- 26.5; k <- function(c, base = 0) (base + slope * c) * (1 + 0.02 * (T0 - 25)) + rnorm(length(c), 0, 1.2)
    st <- data.frame(conc_mg_l = c(0, 50, 100, 200, 300, 400, 500), temp_c = T0); st$cond_us_cm <- round(k(st$conc_mg_l, 3), 1)
    bl <- data.frame(cond_us_cm = round(k(rep(0, 8), 3), 1)); un <- data.frame(cond_us_cm = round(k(rep(240, 3), 3), 1), temp_c = T0)
    ad <- data.frame(added_mg_l = c(0, 50, 100, 150, 200), temp_c = T0); ad$cond_us_cm <- round(k(120 + ad$added_mg_l, 3) * 0.93, 1)      # matrix suppresses the response by 7 percent
    cat("Practice data: 7 salt standards (0 to 500 mg/L, true slope 2.0 uS/cm per mg/L), 8 blanks, an unknown of 240 mg/L, and a standard addition on a sample of 120 mg/L\n\n") }
  st$k25 <- to25(st$cond_us_cm, st$temp_c); fit <- lm(k25 ~ conc_mg_l, st); b <- coef(fit)
  cat(sprintf("Calibration: conductivity(25 C) = %.2f + %.3f x concentration; R squared %.4f; residual sd %.2f uS/cm\n", b[1], b[2], summary(fit)$r.squared, sigma(fit)))
  sb <- sd(to25(bl$cond_us_cm, 25)); lod <- 3 * sb / b[2]; loq <- 10 * sb / b[2]
  cat(sprintf("Blank sd %.2f uS/cm (n = %d): detection limit (3 s / slope) = %.1f mg/L; quantitation limit (10 s / slope) = %.1f mg/L\n", sb, nrow(bl), lod, loq))
  y0 <- mean(to25(un$cond_us_cm, un$temp_c)); x0 <- (y0 - b[1]) / b[2]
  se <- sigma(fit) * sqrt(1 / nrow(un) + 1 / nrow(st) + (x0 - mean(st$conc_mg_l))^2 / sum((st$conc_mg_l - mean(st$conc_mg_l))^2)) / b[2]
  cat(sprintf("Unknown: %.1f mg/L with an approximate 95%% interval of %.1f to %.1f mg/L (n = %d)\n", x0, x0 - qt(.975, nrow(st) - 2) * se, x0 + qt(.975, nrow(st) - 2) * se, nrow(un)))
  ad$k25 <- to25(ad$cond_us_cm, ad$temp_c); fa <- lm(k25 ~ added_mg_l, ad); ba <- coef(fa)
  cat(sprintf("\nStandard addition: slope %.3f (external calibration slope %.3f, ratio %.2f); the line meets the axis at %.1f mg/L in the sample\n", ba[2], b[2], ba[2] / b[2], ba[1] / ba[2]))
  cat(sprintf("Reading the same sample off the external calibration instead gives %.1f mg/L\n", (ba[1] - b[1]) / b[2]))
}
