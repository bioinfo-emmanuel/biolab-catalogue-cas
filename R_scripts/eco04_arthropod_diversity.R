# ECO-04  Leaf-litter arthropod extraction: richness, diversity indices, and comparing sites. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript eco04_arthropod_diversity.R [catch.csv]
# catch.csv: site, sample, morphospecies, count   (one row per morphospecies per litter sample)
set.seed(4)
shannon <- function(x) { p <- x[x > 0] / sum(x); -sum(p * log(p)) }
simpson <- function(x) { p <- x / sum(x); 1 - sum(p^2) }
pielou <- function(x) shannon(x) / log(sum(x > 0))
bray_curtis <- function(a, b) sum(abs(a - b)) / sum(a + b)
rarefy <- function(x, m) { x <- x[x > 0]; N <- sum(x); if (m > N) NA else sum(1 - exp(lchoose(N - x, m) - lchoose(N, m))) }

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) d <- read.csv(args[1]) else {
    sp <- paste0("morpho", 1:25); mk <- function(site, ev, n_s) do.call(rbind, lapply(1:n_s, function(s) { p <- ev^(0:24); p <- p / sum(p); cnt <- rmultinom(1, rpois(1, 60), p)[, 1]; data.frame(site = site, sample = s, morphospecies = sp[cnt > 0], count = cnt[cnt > 0]) }))
    d <- rbind(mk("forest_litter", 0.90, 6), mk("garden_litter", 0.72, 6)); cat("Practice data: 6 litter samples from each of two invented sites (about 60 animals per sample; 25 morphospecies in the pool)\n\n") }
  tab <- tapply(d$count, list(paste(d$site, d$sample), d$morphospecies), sum, default = 0); site_of <- sub(" .*", "", rownames(tab))
  per <- data.frame(site = site_of, animals = rowSums(tab), richness = rowSums(tab > 0), H = apply(tab, 1, shannon), D = apply(tab, 1, simpson), J = apply(tab, 1, pielou))
  cat("Per-sample means and sd (samples are the replicates):\n"); ag <- do.call(data.frame, aggregate(cbind(animals, richness, H, D, J) ~ site, per, function(x) c(mean = mean(x), sd = sd(x)))); ag[, -1] <- round(ag[, -1], 2); print(ag, row.names = FALSE)
  pool <- rowsum(tab, site_of); cat("\nPooled by site:\n"); print(data.frame(site = rownames(pool), animals = rowSums(pool), richness = rowSums(pool > 0), H = round(apply(pool, 1, shannon), 2), D = round(apply(pool, 1, simpson), 3), J = round(apply(pool, 1, pielou), 2)), row.names = FALSE)
  m <- min(rowSums(pool)); cat(sprintf("\nRarefied richness at %d animals (the smaller pooled catch): ", m)); print(round(apply(pool, 1, rarefy, m = m), 1))
  s1 <- unique(site_of)[1]; s2 <- unique(site_of)[2]; cat(sprintf("Bray-Curtis dissimilarity of the pooled catches: %.2f (0 = identical, 1 = no shared morphospecies)\n", bray_curtis(pool[s1, ], pool[s2, ])))
  for (v in c("richness", "H")) { w <- wilcox.test(per[[v]][per$site == s1], per[[v]][per$site == s2], exact = FALSE); cat(sprintf("Samples as replicates, %s: Wilcoxon p = %.4f\n", v, w$p.value)) }
  cat("\nSample-size effect on richness: mean number of morphospecies found with 1 to 6 pooled samples (site 1)\n"); sub <- tab[site_of == s1, ]; for (k in 1:6) cat(sprintf("%d sample(s): %.1f\n", k, mean(replicate(200, sum(colSums(sub[sample(nrow(sub), k), , drop = FALSE]) > 0)))))
}
