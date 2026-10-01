# ECO-01  Quadrat sampling: density, frequency, cover, and sample size. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript eco01_quadrat_sampling.R [quadrats.csv]
# quadrats.csv: habitat, quadrat, species, count, cover_pct   (one row per species per quadrat; quadrat area is set below; include empty quadrats as species "none", count 0)
set.seed(11)
quadrat_area_m2 <- 1

summ <- function(d, area = quadrat_area_m2) {
  nq <- length(unique(d$quadrat)); sp <- split(d[d$species != "none", ], d[d$species != "none", ]$species)
  r <- do.call(rbind, lapply(names(sp), function(s) { x <- sp[[s]]; data.frame(species = s, density_per_m2 = sum(x$count) / (nq * area), frequency_pct = 100 * length(unique(x$quadrat[x$count > 0])) / nq, mean_cover_pct = sum(x$cover_pct) / nq) }))
  r$rel_density <- 100 * r$density_per_m2 / sum(r$density_per_m2); r$rel_frequency <- 100 * r$frequency_pct / sum(r$frequency_pct); r$rel_cover <- 100 * r$mean_cover_pct / sum(r$mean_cover_pct)
  r$importance_value <- r$rel_density + r$rel_frequency + r$rel_cover; r[order(-r$importance_value), ] }
running_cv <- function(x) sapply(seq(3, length(x)), function(k) sd(x[1:k]) / mean(x[1:k]) / sqrt(k))     # standard error of the mean as a fraction of the mean

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) d <- read.csv(args[1]) else {
    mk <- function(h, n, lam) do.call(rbind, lapply(1:n, function(q) do.call(rbind, lapply(names(lam), function(s) { c <- rpois(1, lam[[s]]); data.frame(habitat = h, quadrat = q, species = s, count = c, cover_pct = ifelse(c == 0, 0, round(min(100, c * runif(1, 0.8, 2.4)), 0))) }))))
    d <- rbind(mk("grassland", 30, c(grass_A = 14, herb_B = 3, herb_C = 1.2, sedge_D = 0.4)), mk("shrub_edge", 30, c(grass_A = 4, herb_B = 1.5, herb_C = 2.8, shrub_E = 1.0)))
    cat("Practice data: 30 quadrats of 1 m2 in each of two invented habitats\n\n") }
  for (h in unique(d$habitat)) { cat("Habitat:", h, "\n"); print(format(summ(d[d$habitat == h, ]), digits = 3), row.names = FALSE); cat("\n") }
  tot <- aggregate(count ~ habitat + quadrat, d, sum); mu <- tapply(tot$count, tot$habitat, mean); sdv <- tapply(tot$count, tot$habitat, sd)
  for (h in names(mu)) cat(sprintf("Total plants per m2 in %s: mean %.1f (sd %.1f)\n", h, mu[[h]], sdv[[h]]))
  w <- wilcox.test(count ~ habitat, tot, exact = FALSE); t <- t.test(count ~ habitat, tot); cat(sprintf("Welch t test p = %.2e; Wilcoxon p = %.2e\n", t$p.value, w$p.value))
  cat("\nHow many quadrats? Standard error of the mean density as a fraction of the mean, by number of quadrats (first habitat):\n")
  x <- tot$count[tot$habitat == unique(tot$habitat)[1]]; cv <- running_cv(x); for (k in c(5, 10, 15, 20, 30)) if (k <= length(x)) cat(sprintf("n = %2d  se/mean = %.2f\n", k, cv[k - 2]))
  cat(sprintf("Quadrats needed for se/mean = 0.10 (from the pilot's mean and sd): %.0f\n", ceiling((sd(x) / mean(x) / 0.10)^2)))
  cat("\nQuadrat size: the same population sampled with 0.25 m2 quadrats shows a larger coefficient of variation between quadrats than with 1 m2 quadrats.\n")
  lam <- 20; for (a in c(0.25, 1, 4)) { z <- rpois(2000, lam * a); cat(sprintf("area %.2f m2: mean %.1f, cv %.2f\n", a, mean(z), sd(z) / mean(z))) }
}
