# BOT-01  Stomatal density and index from leaf impressions. Base R only: no packages to install.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
#
# YOUR DATA: a CSV with one row per microscope field and these columns
#   species, leaf_id, surface (upper or lower), field_id, stomata, epidermal, fov_diameter_mm
#   stomata          stomata counted in the field
#   epidermal        ordinary epidermal cells counted in the field (not counting stomata)
#   fov_diameter_mm  diameter of the circular field of view at the objective used, in mm
# The leaf, not the field, is the replicate.
#
# Run in RStudio: set data_file below and click Source. From a terminal: Rscript stomata_analysis.R mydata.csv
# With no file the script makes PRACTICE data. Practice data are simulated. Never report them as results.

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1) args[1] else NA

if (is.na(data_file)) {
  message("No data file given: using simulated PRACTICE data.")
  set.seed(21)
  species <- c("Sp1", "Sp2", "Sp3", "Sp4")
  lower_mean <- c(Sp1 = 300, Sp2 = 200, Sp3 = 120, Sp4 = 250)   # stomata per mm2, lower surface
  upper_mean <- c(Sp1 = 30, Sp2 = 90, Sp3 = 0, Sp4 = 140)
  rows <- list(); fov <- 0.45
  area <- pi * (fov / 2)^2
  for (sp in species) for (leaf in 1:6) for (surf in c("upper", "lower")) for (fld in 1:5) {
    lam <- (if (surf == "lower") lower_mean[sp] else upper_mean[sp]) * area * rlnorm(1, 0, 0.05)
    leaf_eff <- rnorm(1, 1, 0.08)
    s <- rpois(1, max(lam * leaf_eff, 0)); e <- rpois(1, 900 * area)
    rows[[length(rows) + 1]] <- data.frame(species = sp, leaf_id = paste0(sp, "_", leaf), surface = surf,
                                            field_id = fld, stomata = s, epidermal = e, fov_diameter_mm = fov)
  }
  d <- do.call(rbind, rows)
} else {
  d <- read.csv(data_file, stringsAsFactors = FALSE)
}

d$area_mm2 <- pi * (d$fov_diameter_mm / 2)^2
d$density <- d$stomata / d$area_mm2                          # stomata per mm2
d$index <- 100 * d$stomata / (d$stomata + d$epidermal)       # stomatal index, percent

# The leaf is the replicate: average the fields within each leaf and surface first.
leaf <- aggregate(cbind(density, index) ~ species + leaf_id + surface, data = d, FUN = mean)

cat("\n1. Mean stomatal density (per mm2) and index (%) by species and surface, from leaf means\n")
print(aggregate(cbind(density, index) ~ species + surface, data = leaf, FUN = function(x) round(mean(x), 1)), row.names = FALSE)

cat("\n2. Upper versus lower surface, paired by leaf (Holm-adjusted p across species)\n")
res <- lapply(split(leaf, leaf$species), function(x) {
  w <- reshape(x[, c("leaf_id", "surface", "density")], idvar = "leaf_id", timevar = "surface", direction = "wide")
  tt <- t.test(w$density.lower, w$density.upper, paired = TRUE)
  data.frame(mean_lower_minus_upper = unname(tt$estimate), t = unname(tt$statistic), df = unname(tt$parameter), p = tt$p.value)
})
tab <- do.call(rbind, res); tab$species <- rownames(tab); tab$p_holm <- p.adjust(tab$p, "holm")
tab[, c("mean_lower_minus_upper", "t", "p", "p_holm")] <- round(tab[, c("mean_lower_minus_upper", "t", "p", "p_holm")], 4)
print(tab[, c("species", "mean_lower_minus_upper", "t", "df", "p", "p_holm")], row.names = FALSE)

cat("\n3. Differences among species, lower surface density (one-way ANOVA on leaf means)\n")
lw <- leaf[leaf$surface == "lower", ]
fit <- aov(density ~ species, data = lw)
print(summary(fit))
print(TukeyHSD(fit))

par(mfrow = c(1, 2))
boxplot(density ~ species + surface, data = leaf, las = 2, ylab = "Stomata per mm2", main = "Density", cex.axis = 0.7)
boxplot(index ~ species + surface, data = leaf, las = 2, ylab = "Stomatal index (%)", main = "Index", cex.axis = 0.7)
