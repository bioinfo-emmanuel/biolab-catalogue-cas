# AP-05  Hollow organ wall structure: layer thickness ratios from prepared slides. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript ap05_wall_layers.R [layers.csv]
# layers.csv: organ, field, layer, thickness   (thickness in any single unit, measured along a line perpendicular to the wall, for example in mm on a photograph)
# Ratios and proportions of wall thickness need no calibration.
set.seed(5)
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) d <- read.csv(args[1]) else {
    spec <- list(artery = c(intima = 1, media = 6, adventitia = 2.5), vein = c(intima = 0.8, media = 1.6, adventitia = 2.2), small_intestine = c(mucosa = 4, submucosa = 1.5, muscularis = 2, serosa = 0.4))
    d <- do.call(rbind, lapply(names(spec), function(o) do.call(rbind, lapply(1:10, function(f) data.frame(organ = o, field = f, layer = names(spec[[o]]), thickness = round(spec[[o]] * rlnorm(length(spec[[o]]), 0, 0.15) * 10, 1))))))
    cat("Practice data: 10 fields for each of three invented organs (units invented); layer proportions are made up for practice, not measurements of real tissues.\n\n") }
  tot <- aggregate(thickness ~ organ + field, d, sum); names(tot)[3] <- "wall"; d <- merge(d, tot); d$prop <- d$thickness / d$wall
  res <- aggregate(prop ~ organ + layer, d, function(x) c(mean = mean(x), sd = sd(x))); res <- do.call(data.frame, res); names(res)[3:4] <- c("mean_prop", "sd_prop"); res$mean_pct <- round(100 * res$mean_prop, 1); print(res[order(res$organ, -res$mean_prop), c("organ", "layer", "mean_pct", "sd_prop")], digits = 3, row.names = FALSE)
  cat("\nWall thickness (units) by organ, mean and sd over fields:\n"); ag <- do.call(data.frame, aggregate(wall ~ organ, tot, function(x) c(mean = mean(x), sd = sd(x)))); ag[, -1] <- round(ag[, -1], 1); print(ag, row.names = FALSE)
  m <- d[d$layer == "media" & d$organ %in% c("artery", "vein"), ]; if (nrow(m)) { w <- wilcox.test(prop ~ organ, m, exact = FALSE); cat(sprintf("\nProportion of wall that is media, artery vs vein: Wilcoxon p = %.2e\n", w$p.value)) }
  cat("\nHow many fields are enough? sd of the mean media proportion for the artery by number of fields (bootstrap)\n"); x <- d$prop[d$organ == "artery" & d$layer == "media"]; for (n in c(3, 5, 10)) cat(sprintf("n = %2d  se = %.4f\n", n, sd(replicate(500, mean(sample(x, n, replace = TRUE))))))
  cat("Section thickness and tilt: a wall cut obliquely looks thicker, so measure only where the layers run straight and the lumen is round.\n")
}
