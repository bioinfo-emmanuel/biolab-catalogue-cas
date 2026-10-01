# CVA-04  Gill respiratory surface area and body size in fish. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript cva04_gill_area.R [fish.csv]
# fish.csv: fish, mass_g, total_filament_mm (all filaments on all gill arches, both sides of the head), lamellae_per_mm (one side of a filament),
#           lamella_length_mm, lamella_height_mm  (mean of about 10 lamellae spread along the filaments)
# Formula (Hughes 1984, J Mar Biol Assoc UK 64:637-655): total lamellae = 2 x L x n; gill area A = total lamellae x bilateral area of an average lamella.
# The bilateral area is 2 x (0.5 x length x height) when the lamella is treated as a triangle-like plate; adjust shape_factor to your own outline tracing.
set.seed(88)
gill_area <- function(L_mm, n_per_mm, bil_area_mm2) 2 * L_mm * n_per_mm * bil_area_mm2
bilateral_area <- function(len, ht, shape_factor = 0.5) 2 * shape_factor * len * ht

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  cat("Check against the icefish measurements in Hughes (1972) J Exp Biol 56(2):481, Table 1 (Chaenocephalus aceratus): L = 15,904 mm, 9.7 lamellae per mm on one side, bilateral area 0.300 mm2\n")
  cat(sprintf("total lamellae = %.0f (published 308,537); gill area = %.0f mm2 (published 92,561)\n\n", 2 * 15904 * 9.7, gill_area(15904, 9.7, 0.300)))
  if (length(args) >= 1) d <- read.csv(args[1]) else {
    n <- 12; m <- round(exp(runif(n, log(80), log(600))), 0)
    d <- data.frame(fish = 1:n, mass_g = m, total_filament_mm = round(230 * m^0.55 * rnorm(n, 1, 0.05), 0), lamellae_per_mm = round(rnorm(n, 20, 0.8) * (m / 300)^-0.05, 1),
                    lamella_length_mm = round(0.28 * (m / 300)^0.30 * rnorm(n, 1, 0.04), 3), lamella_height_mm = round(0.10 * (m / 300)^0.10 * rnorm(n, 1, 0.04), 3))
    cat("Practice data: 12 simulated fish from 80 to 600 g. The dimensions are invented for practice, not measurements of any species.\n\n") }
  d$bil_area <- bilateral_area(d$lamella_length_mm, d$lamella_height_mm); d$total_lamellae <- 2 * d$total_filament_mm * d$lamellae_per_mm
  d$gill_area_mm2 <- gill_area(d$total_filament_mm, d$lamellae_per_mm, d$bil_area); d$area_per_g <- d$gill_area_mm2 / d$mass_g
  print(round(d[, c("fish", "mass_g", "total_filament_mm", "total_lamellae", "bil_area", "gill_area_mm2", "area_per_g")], 3), row.names = FALSE)
  cat("\nAllometry: log10(gill area) against log10(body mass)\n"); f <- lm(log10(gill_area_mm2) ~ log10(mass_g), d); ci <- confint(f)[2, ]
  cat(sprintf("slope %.3f (95%% CI %.3f to %.3f), R squared %.3f. A slope of 1 would mean gill area rises in proportion to mass; less than 1 means area per gram falls as fish grow.\n", coef(f)[2], ci[1], ci[2], summary(f)$r.squared))
  for (v in c("total_filament_mm", "lamellae_per_mm", "bil_area")) { g <- lm(log10(d[[v]]) ~ log10(d$mass_g)); cat(sprintf("component %-18s slope %.3f\n", v, coef(g)[2])) }
  cat("\nSensitivity: a 10 percent error in lamellar length changes bilateral area and hence gill area by 10 percent; the same error in lamellae per mm changes area by 10 percent.\n")
}
