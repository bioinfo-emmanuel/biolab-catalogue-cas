# ECO-03  Leaf-litter decomposition rates with litter bags. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript eco03_litter_decay.R [bags.csv]
# bags.csv: species, site, bag, week, m0_dry_g (estimated starting dry mass of that bag), mt_dry_g (dry mass of the litter recovered)
# Starting dry mass: weigh the air-dry litter, and oven-dry a separate subsample to get the ratio of dry to air-dry mass.
set.seed(3)
decay_k <- function(m0, mt, weeks) -log(mt / m0) / weeks                       # per week, from M(t) = M0 exp(-k t)
fit_k <- function(d) { f <- lm(log(mt_dry_g / m0_dry_g) ~ 0 + week, d); k <- -coef(f)[[1]]; ci <- -rev(as.numeric(confint(f)[1, ])); c(k = k, lower = ci[1], upper = ci[2], half_life_weeks = log(2) / k) }

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) d <- read.csv(args[1]) else {
    k <- c(species_1 = 0.11, species_2 = 0.04); weeks <- c(2, 4, 6, 8); bags <- 4
    d <- do.call(rbind, lapply(names(k), function(sp) do.call(rbind, lapply(weeks, function(w) { m0 <- rnorm(bags, 10, 0.4); data.frame(species = sp, site = "site_A", bag = paste0(sp, "_w", w, "_", 1:bags), week = w, m0_dry_g = round(m0, 1), mt_dry_g = round(m0 * exp(-k[[sp]] * w * rnorm(bags, 1, 0.12)), 1)) }))))
    cat("Practice data: two invented litter types, 4 bags at each of 4 retrieval times (weeks 2, 4, 6, 8), 10 g dry at the start, balance reading to 0.1 g\n\n") }
  d$prop_remaining <- d$mt_dry_g / d$m0_dry_g; d$k_bag <- decay_k(d$m0_dry_g, d$mt_dry_g, d$week)
  cat("Mass remaining (percent of the starting mass) by retrieval week:\n"); print(round(100 * tapply(d$prop_remaining, list(d$species, d$week), mean), 1))
  res <- do.call(rbind, lapply(split(d, d$species), function(x) { r <- fit_k(x); data.frame(species = x$species[1], k_per_week = r["k"], ci_low = r["lower"], ci_high = r["upper"], half_life_weeks = r["half_life_weeks"], bags = nrow(x)) })); print(format(res, digits = 3), row.names = FALSE)
  m <- lm(log(prop_remaining) ~ 0 + week:species, d); m0 <- lm(log(prop_remaining) ~ 0 + week, d); cat(sprintf("\nDo the two litter types decay at different rates? F test of separate slopes against one common slope: p = %.2e\n", anova(m0, m)$`Pr(>F)`[2]))
  cat("Check: 10.0 g falling to 6.7 g in 4 weeks gives k =", round(decay_k(10, 6.7, 4), 3), "per week, half-life", round(log(2) / decay_k(10, 6.7, 4), 1), "weeks\n")
  cat("Reading error: a 5 g bag read to 0.1 g has about", round(100 * 0.1 / sqrt(12) / 5, 1), "percent reading error, and the recovered mass of a heavily decomposed bag (1 g) has", round(100 * 0.1 / sqrt(12) / 1, 1), "percent.\n")
}
