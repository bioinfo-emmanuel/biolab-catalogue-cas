# BCH-04  Salting out proteins: ammonium sulfate amounts and fractional precipitation of egg white. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript bch04_ammonium_sulfate.R [pellets.csv]
# pellets.csv: fraction, dry_pellet_g   (dry mass of heat-coagulated protein recovered from each redissolved pellet; the raw pellet also holds salt and water, so it is not weighed as protein)
# The gram amounts use the equation of Wood (1976)-type tables in the form used by EnCor Biotechnology:
#   G = Sat x (S2 - S1) / (1 - Vs x Sat x S2 / 1000)   grams per litre, with Sat the grams that saturate 1 L of solution
# and Vs the apparent specific volume of the salt (mL per g). Constants at 20 C from Wingfield (Curr Protoc Protein Sci table A.3F.1): Sat 536.34 g/L, Vs 0.5414 mL/g.
Sat <- 536.34; Vs <- 0.5414
grams_to_add <- function(S1, S2, volume_ml, Sat = 536.34, Vs = 0.5414) volume_ml / 1000 * Sat * (S2 - S1) / (1 - Vs * Sat * S2 / 1000)
set.seed(44)

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  cat("Grams of solid ammonium sulfate to add per 100 mL of solution at 20 C (S1 -> S2 as fractions of saturation)\n")
  tab <- expand.grid(S1 = c(0, 0.2, 0.4), S2 = c(0.4, 0.5, 0.6, 0.8, 1.0)); tab <- tab[tab$S2 > tab$S1, ]; tab$g_per_100ml <- round(grams_to_add(tab$S1, tab$S2, 100), 1); print(tab, row.names = FALSE)
  cat(sprintf("\nCheck against the 0 C table of Wood (1976) as printed in Cold Spring Harbor Protocols: 0 to 50%% at 0 C is 29.5 g per 100 mL (0 C constants: Sat 514.72, Vs 0.5262); this script gives %.1f\n", grams_to_add(0, 0.5, 100, 514.72, 0.5262)))
  cat(sprintf("Check: 0 to 100%% at 0 C is 70.7 g per 100 mL in the same table; this script gives %.1f\n", grams_to_add(0, 1, 100, 514.72, 0.5262)))
  cat("\nAmounts to weigh for one group: 50 mL egg white plus 200 mL water (250 mL), cut at 0 to 33%, then 33 to 50%, then 50 to 100% saturation\n")
  vol <- 250; cuts <- c(0, 0.33, 0.5, 1.0); tot <- 0
  for (i in 1:3) { g <- grams_to_add(cuts[i], cuts[i + 1], vol); vol <- vol + g * Vs; tot <- tot + g
    cat(sprintf("Cut %d (%.0f to %.0f%%): add %.1f g; volume after dissolving about %.0f mL (the pellet removed between cuts is small, so the volume is not corrected for it)\n", i, 100 * cuts[i], 100 * cuts[i + 1], g, vol)) }
  cat(sprintf("Total salt per group about %.0f g; for 4 groups about %.0f g\n", tot, 4 * tot))
  d <- if (length(args) >= 1) read.csv(args[1]) else { cat("\nPractice data: masses (g) of heat-coagulated protein recovered from the redissolved pellets of three cuts of 250 mL diluted egg white (50 mL plus 200 mL water), simulated for four groups; the values are invented and are not real egg-white fractions\n"); data.frame(group = rep(1:4, each = 3), fraction = rep(c("0-33%", "33-50%", "50-100%"), 4), dry_pellet_g = round(rep(c(0.9, 0.5, 4.6), 4) * runif(12, 0.85, 1.15), 2)) }
  agg <- aggregate(dry_pellet_g ~ fraction, d, function(x) c(mean = mean(x), sd = sd(x))); agg <- do.call(data.frame, agg); names(agg) <- c("fraction", "mean_g", "sd_g"); agg$share_pct <- 100 * agg$mean_g / sum(agg$mean_g); print(format(agg, digits = 3), row.names = FALSE)
  cat("\nA balance reading to 0.1 g cannot see a fraction of 0.05 g; use large starting volumes and pre-weighed filter paper for the coagulum.\n")
}
