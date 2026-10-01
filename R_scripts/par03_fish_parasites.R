# PAR-03  Parasite survey in fish: prevalence, intensity, abundance, and aggregation. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript par03_fish_parasites.R [fish.csv]
# fish.csv: fish, group (for example site or species), length_cm, parasites (count of one parasite group in that fish)
# Definitions follow Bush AO, Lafferty KD, Lotz JM, Shostak AW (1997) J Parasitol 83:575-583:
#   prevalence = infected hosts / hosts examined; mean intensity = parasites / infected hosts; mean abundance = parasites / all hosts examined.
set.seed(3)
nb_k <- function(x) { nll <- function(p) -sum(dnbinom(x, size = exp(p[2]), mu = exp(p[1]), log = TRUE)); o <- optim(c(log(mean(x) + 0.1), 0), nll, method = "BFGS"); c(mu = exp(o$par[1]), k = exp(o$par[2])) }
indices <- function(x) { inf <- x > 0; c(hosts = length(x), infected = sum(inf), prevalence = mean(inf), mean_intensity = if (any(inf)) mean(x[inf]) else NA, mean_abundance = mean(x), variance_to_mean = var(x) / mean(x), top20_share = { y <- sort(x, decreasing = TRUE); sum(y[seq_len(ceiling(0.2 * length(y)))]) / sum(y) }) }
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) d <- read.csv(args[1]) else {
    mk <- function(g, n, p, mu, k, mean_len) { len <- round(rnorm(n, mean_len, 3), 1); cnt <- ifelse(runif(n) < p, 1 + rnbinom(n, size = k, mu = mu), 0); data.frame(fish = paste0(g, 1:n), group = g, length_cm = len, parasites = cnt) }
    d <- rbind(mk("pond_A", 40, 0.35, 4, 0.6, 15), mk("market_B", 40, 0.62, 7, 0.5, 18)); cat("Practice data: 40 invented fish from each of two sources; the infection levels are invented.\n\n") }
  res <- do.call(rbind, lapply(split(d, d$group), function(x) { r <- indices(x$parasites); ci <- binom.test(r[["infected"]], r[["hosts"]])$conf.int; data.frame(group = x$group[1], t(round(r, 2)), prev_lo = ci[1], prev_hi = ci[2]) })); print(res, digits = 3, row.names = FALSE)
  g <- unique(d$group); tab <- table(d$group, d$parasites > 0); f <- fisher.test(tab); cat(sprintf("\nPrevalence, %s vs %s: Fisher exact p = %.4f, odds ratio %.1f\n", g[1], g[2], f$p.value, f$estimate))
  a <- d[d$group == g[1] & d$parasites > 0, "parasites"]; b <- d[d$group == g[2] & d$parasites > 0, "parasites"]; cat(sprintf("Intensity among infected fish: Wilcoxon p = %.4f (medians %.0f and %.0f)\n", wilcox.test(a, b, exact = FALSE)$p.value, median(a), median(b)))
  cat("\nAggregation (negative binomial k over all fish; small k means strong aggregation):\n"); for (gr in g) { r <- nb_k(d$parasites[d$group == gr]); cat(sprintf("%s: mean %.2f, k = %.2f\n", gr, r["mu"], r["k"])) }
  ct <- cor.test(d$length_cm, d$parasites, method = "spearman", exact = FALSE); cat(sprintf("\nParasites against host length across all fish: Spearman rho = %.2f, p = %.3g\n", ct$estimate, ct$p.value))
  cat("A parasite count is a property of a fish, so the fish is the unit; do not test parasites as if they were independent animals.\n")
  cat("Check: 20 infected of 40 fish is a prevalence of", 100 * 20 / 40, "percent; 120 parasites among 20 infected fish gives a mean intensity of", 120 / 20, "and a mean abundance of", 120 / 40, "\n")
}
