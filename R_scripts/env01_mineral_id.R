# ENV-01  Mineral identification: hardness, streak, density, and acid test. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript env01_mineral_id.R [observations.csv]
# observations.csv: sample, hardness (Mohs, your estimate), streak, mass_air_g, mass_water_g (mass when suspended in water, or NA),
#                   volume_ml (by displacement, or NA), fizzes_acid (yes/no/powder_only)
# The reference table holds typical values for common minerals. Check it against your own textbook or database before publishing.
set.seed(55)
ref <- data.frame(mineral = c("quartz", "calcite", "orthoclase feldspar", "gypsum", "halite", "muscovite", "biotite", "galena", "pyrite", "hematite", "magnetite", "fluorite", "talc"),
  h_lo = c(7, 3, 6, 2, 2.5, 2, 2.5, 2.5, 6, 5, 5.5, 4, 1), h_hi = c(7, 3, 6, 2, 2.5, 2.5, 3, 2.5, 6.5, 6, 6.5, 4, 1),
  sg_lo = c(2.65, 2.71, 2.55, 2.30, 2.17, 2.76, 2.7, 7.4, 4.9, 4.9, 5.15, 3.0, 2.7), sg_hi = c(2.65, 2.71, 2.63, 2.33, 2.17, 3.0, 3.4, 7.6, 5.2, 5.3, 5.2, 3.3, 2.8),
  streak = c("white", "white", "white", "white", "white", "white", "white", "gray", "green-black", "red-brown", "black", "white", "white"),
  acid = c("no", "yes", "no", "no", "no", "no", "no", "no", "no", "no", "no", "no", "no"), stringsAsFactors = FALSE)

density_from_masses <- function(m_air, m_water_susp = NA, volume_ml = NA) ifelse(!is.na(m_water_susp), m_air / (m_air - m_water_susp), ifelse(!is.na(volume_ml), m_air / volume_ml, NA))
score <- function(obs, ref) {
  s <- sapply(seq_len(nrow(ref)), function(i) { r <- ref[i, ]; pts <- 0
    if (!is.na(obs$hardness)) pts <- pts + (obs$hardness >= r$h_lo - 0.6 & obs$hardness <= r$h_hi + 0.6)
    if (!is.na(obs$sg)) pts <- pts + 2 * (obs$sg >= r$sg_lo * 0.93 & obs$sg <= r$sg_hi * 1.07)
    if (!is.na(obs$streak)) pts <- pts + (tolower(obs$streak) == r$streak)
    if (!is.na(obs$fizzes_acid) && obs$fizzes_acid %in% c("yes", "no")) pts <- pts + (obs$fizzes_acid == r$acid)
    pts })
  data.frame(mineral = ref$mineral, score = s)[order(-s), ] }

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) obs <- read.csv(args[1], stringsAsFactors = FALSE) else {
    obs <- data.frame(sample = c("A", "B", "C", "D"), hardness = c(7, 3, 2.5, 6), streak = c("white", "white", "gray", "green-black"),
      mass_air_g = c(38.2, 27.9, 61.4, 44.0), mass_water_g = c(23.4, 17.6, 53.2, 35.3), volume_ml = NA, fizzes_acid = c("no", "yes", "no", "no"), stringsAsFactors = FALSE)
    cat("Practice data: four invented specimens (quartz-like, calcite-like, galena-like and pyrite-like properties)\n\n") }
  obs$sg <- density_from_masses(obs$mass_air_g, obs$mass_water_g, obs$volume_ml); print(obs[, c("sample", "hardness", "streak", "sg", "fizzes_acid")], digits = 3, row.names = FALSE)
  cat("\nBest matches:\n"); for (i in seq_len(nrow(obs))) { sc <- score(obs[i, ], ref); cat(sprintf("%s: %s (score %d), then %s (%d)\n", obs$sample[i], sc$mineral[1], sc$score[1], sc$mineral[2], sc$score[2])) }
  cat("\nHow much does a balance reading to 0.1 g matter for specific gravity? A 5 g sample of quartz (true 2.65): m_air 5.0, suspended 3.1 gives", round(density_from_masses(5.0, 3.1), 2), "; with each reading off by 0.1 g the value ranges over",
      paste(round(range(density_from_masses(5.0 + c(-.1, .1), 3.1 + c(.1, -.1))), 2), collapse = " to "), "\n")
  cat("A 50 g sample with 0.1 g reading error on each mass changes specific gravity by roughly", round(100 * 0.2 / (50 - 31.13), 1), "percent for quartz (mass in water 31.1 g).\n")
}
