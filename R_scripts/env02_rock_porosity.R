# ENV-02  Porosity and density of rocks by water absorption. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript env02_rock_porosity.R [rocks.csv]
# rocks.csv: rock, replicate, m_dry_g, m_sat_g (saturated surface-dry mass), volume_ml (bulk volume by water displacement of the saturated sample)
# Formulas follow the usual saturation and calliper-style method: absorption = (m_sat - m_dry) / m_dry; bulk density = m_dry / V;
# apparent porosity = (m_sat - m_dry) / (rho_water x V). Set water density for your temperature (about 0.997 g/mL at 25 C).
set.seed(66)
rho_w <- 0.997
props <- function(m_dry, m_sat, v) data.frame(absorption_pct = 100 * (m_sat - m_dry) / m_dry, bulk_density = m_dry / v, porosity_pct = 100 * (m_sat - m_dry) / (rho_w * v))

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) d <- read.csv(args[1]) else {
    types <- data.frame(rock = c("limestone", "sandstone", "andesite", "basalt", "coral limestone"), rho = c(2.55, 2.20, 2.55, 2.85, 1.85), phi = c(0.05, 0.18, 0.05, 0.03, 0.30))
    d <- do.call(rbind, lapply(seq_len(nrow(types)), function(i) { v <- rnorm(3, 100, 6); md <- types$rho[i] * v * rnorm(3, 1, 0.01); data.frame(rock = types$rock[i], replicate = 1:3, m_dry_g = round(md, 1), m_sat_g = round(md + types$phi[i] * v * rho_w, 1), volume_ml = round(v, 0)) }))
    cat("Practice data: five invented rock types, three pieces of about 100 mL each; masses read to 0.1 g. The porosities and densities are made up for practice, not measured values for these rocks.\n\n") }
  p <- cbind(d, props(d$m_dry_g, d$m_sat_g, d$volume_ml))
  res <- do.call(rbind, lapply(split(p, p$rock), function(x) data.frame(rock = x$rock[1], n = nrow(x), absorption_pct = mean(x$absorption_pct), sd_abs = sd(x$absorption_pct), bulk_density = mean(x$bulk_density), porosity_pct = mean(x$porosity_pct), sd_por = sd(x$porosity_pct))))
  print(format(res, digits = 3), row.names = FALSE)
  cat("\nRelationship between bulk density and porosity across the samples:\n"); f <- lm(porosity_pct ~ bulk_density, p); print(round(summary(f)$coefficients, 3)); cat(sprintf("R squared %.2f\n", summary(f)$r.squared))
  cat("\nReading error: a 100 g piece with 0.1 g reading error on each mass, absorbing 1 g (1 percent): the absorbed mass is uncertain by about", round(100 * sqrt(2) * 0.1 / 1, 0), "percent of itself.\n")
  cat("Low-porosity rocks (below about 1 percent absorption) need larger pieces, or a balance that reads to 0.01 g.\n")
}
