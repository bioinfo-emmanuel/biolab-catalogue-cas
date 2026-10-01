# PAR-06  Overdispersion of parasite counts: Poisson versus negative binomial. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript par06_negative_binomial.R [counts.csv]
# counts.csv column: count (parasites per host). With no argument the script uses R's built-in InsectSprays data
# (counts of insects per plot) as an aggregated count example, and a simulated host-parasite sample.

nb_fit <- function(x) {                         # maximum likelihood for mean mu and aggregation parameter k (size)
  nll <- function(par) -sum(dnbinom(x, size = exp(par[2]), mu = exp(par[1]), log = TRUE))
  o <- optim(c(log(mean(x)), 0), nll, method = "BFGS"); c(mu = exp(o$par[1]), k = exp(o$par[2]), loglik = -o$value)
}
pois_ll <- function(x) sum(dpois(x, mean(x), log = TRUE))
compare <- function(x) { nb <- nb_fit(x); pl <- pois_ll(x); n <- length(x)
  c(mean = mean(x), variance = var(x), mu = nb[["mu"]], k = nb[["k"]], ll_pois = pl, ll_nb = nb[["loglik"]],
    AIC_pois = -2 * pl + 2, AIC_nb = -2 * nb[["loglik"]] + 4, LRT_p = pchisq(2 * (nb[["loglik"]] - pl), 1, lower.tail = FALSE)) }
# Share of all parasites carried by the most heavily infected fraction of hosts.
top_share <- function(x, frac = 0.2) { x <- sort(x, decreasing = TRUE); sum(x[seq_len(ceiling(frac * length(x)))]) / sum(x) }

if (sys.nframe() == 0) {
  set.seed(6)
  args <- commandArgs(trailingOnly = TRUE)
  x <- if (length(args) >= 1) read.csv(args[1])$count else rnbinom(120, size = 0.4, mu = 6)
  if (length(args) < 1) cat("Practice sample: 120 simulated hosts, true mean 6 and true k = 0.4\n")
  cat("\nParasite counts per host\n"); print(round(compare(x), 3))
  cat(sprintf("Share of parasites in the most heavily infected 20%% of hosts: %.0f%% (a Poisson distribution with this mean would give about %.0f%%)\n", 100 * top_share(x), 100 * top_share(rpois(1e5, mean(x)))))
  cat(sprintf("Prevalence (hosts with at least one parasite): %.2f; Poisson would predict %.2f\n", mean(x > 0), 1 - exp(-mean(x))))
  cat("\nBuilt-in InsectSprays data, spray C (12 plots) and spray A (12 plots) counts as a second example\n")
  ins <- datasets::InsectSprays; for (sp in c("A", "C")) { y <- ins$count[ins$spray == sp]; r <- compare(y); cat(sprintf("spray %s: mean %.1f variance %.1f  variance/mean %.2f  k = %.1f\n", sp, r["mean"], r["variance"], r["variance"] / r["mean"], r["k"])) }
  cat("\nHow well does k recover the truth? 200 samples of size 30 and 300 (true k = 0.4)\n")
  for (n in c(30, 300)) { ks <- replicate(200, nb_fit(rnbinom(n, size = 0.4, mu = 6))["k"]); cat(sprintf("n = %3d  median k %.2f  5th to 95th percentile %.2f to %.2f\n", n, median(ks), quantile(ks, .05), quantile(ks, .95))) }
}
