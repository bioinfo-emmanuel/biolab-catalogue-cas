# ACH-01  Gravimetric calibration of glassware. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript ach01_gravimetric_calibration.R [data.csv]
# data.csv columns: device, nominal_ml, mass_g, temp_c   (one row per delivery of water; mass of the water delivered)
# With no argument the script simulates three devices so it runs offline.
set.seed(101)

# Density of air-free water, Tanaka et al. (2001) Metrologia 38:301-309, valid 0 to 40 C, in g/mL.
water_density <- function(t) 0.99997495 * (1 - ((t - 3.983035)^2 * (t + 301.797)) / (522528.9 * (t + 69.34881)))
# Delivered volume in mL from the mass of water in air, with the air-buoyancy correction:
# V = m / (rho_water - rho_air) * (1 - rho_air / rho_weights). rho_air 0.0012 g/mL; balance weights 8.0 g/mL.
delivered_ml <- function(m, t, rho_air = 0.0012, rho_wt = 8.0) m / (water_density(t) - rho_air) * (1 - rho_air / rho_wt)

summarize_device <- function(d, resolution_g = 0.01) {
  v <- delivered_ml(d$mass_g, d$temp_c); nom <- d$nominal_ml[1]; n <- length(v)
  tt <- t.test(v, mu = nom)
  u_bal <- sqrt(2) * resolution_g / sqrt(12) / (water_density(mean(d$temp_c)) - 0.0012)      # tare and gross readings, rectangular distribution, in mL
  c(n = n, nominal = nom, mean = mean(v), sd = sd(v), rsd_pct = 100 * sd(v) / mean(v), bias = mean(v) - nom, bias_pct = 100 * (mean(v) - nom) / nom,
    p_bias = tt$p.value, u_repeat = sd(v) / sqrt(n), u_balance = u_bal) }

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE); res <- 0.01
  if (length(args) >= 1) d <- read.csv(args[1]) else {
    sim <- function(dev, nom, bias, sd, n = 10, temp = 27.5) { v <- rnorm(n, nom + bias, sd); data.frame(device = dev, nominal_ml = nom, mass_g = round(v * (water_density(temp) - 0.0012) / (1 - 0.0012 / 8), 2), temp_c = temp) }
    d <- rbind(sim("cylinder_100mL", 100, -0.8, 0.35), sim("pipette_25mL", 25, 0.02, 0.03), sim("syringe_10mL", 10, 0.25, 0.12))
    cat("Practice data: three simulated devices, 10 deliveries each, water at 27.5 C, balance reading to 0.01 g\n\n") }
  cat("Water density (Tanaka 2001) at 20, 25 and 30 C:", round(water_density(c(20, 25, 30)), 5), "g/mL\n")
  cat("Volume of 1 g of water in air at 20 C (Z factor):", round(delivered_ml(1, 20), 4), "mL/g\n\n")
  out <- do.call(rbind, lapply(split(d, d$device), summarize_device, resolution_g = res)); print(round(out, 4))
  cat("\nEffect of ignoring temperature: reading the volume as mass / 0.9982 (the 20 C density) at 27.5 C\n")
  for (dv in unique(d$device)) { s <- d[d$device == dv, ]; cat(sprintf("%-16s correct mean %.3f mL; with the 20 C density %.3f mL (error %.3f%%)\n", dv, mean(delivered_ml(s$mass_g, s$temp_c)), mean(s$mass_g / 0.9982), 100 * (mean(s$mass_g / 0.9982) / mean(delivered_ml(s$mass_g, s$temp_c)) - 1))) }
  cat("\nReading error of a balance with 0.1 g resolution for a 10 mL delivery (rectangular distribution, tare and gross): ",
      round(100 * sqrt(2) * 0.1 / sqrt(12) / 10, 2), "% of the delivered volume\n", sep = "")
}
