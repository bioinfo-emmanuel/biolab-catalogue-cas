# MIC-05  Antimicrobial resistance surveillance: proportions, intervals, and trends. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript mic05_resistance_surveillance.R [table.csv]
# table.csv columns: organism, drug, year, tested, resistant.  With no argument the script simulates a table
# for two organisms and three drugs over 2015 to 2024 so it runs offline (the numbers are invented).
set.seed(505)

simulate_table <- function() { g <- expand.grid(organism = c("Org_1", "Org_2"), drug = c("Drug_A", "Drug_B", "Drug_C"), year = 2015:2024, stringsAsFactors = FALSE)
  base <- c(Org_1.Drug_A = 0.10, Org_1.Drug_B = 0.30, Org_1.Drug_C = 0.05, Org_2.Drug_A = 0.25, Org_2.Drug_B = 0.20, Org_2.Drug_C = 0.15)
  slope <- c(Org_1.Drug_A = 0.10, Org_1.Drug_B = 0, Org_1.Drug_C = 0.15, Org_2.Drug_A = 0.05, Org_2.Drug_B = -0.10, Org_2.Drug_C = 0)   # log-odds change per year
  k <- paste(g$organism, g$drug, sep = "."); lo <- qlogis(base[k]) + slope[k] * (g$year - 2015); g$tested <- rpois(nrow(g), ifelse(g$organism == "Org_1", 300, 60)); g$resistant <- rbinom(nrow(g), g$tested, plogis(lo)); g }
prop_ci <- function(x, n) { b <- binom.test(x, n); c(prop = x / n, lower = b$conf.int[1], upper = b$conf.int[2]) }
trend_fit <- function(d) { d$yr <- d$year - min(d$year); f <- glm(cbind(resistant, tested - resistant) ~ yr, family = binomial, data = d); s <- summary(f)$coefficients["yr", ]
  ci <- suppressMessages(confint.default(f))["yr", ]; c(OR_per_year = unname(exp(s[1])), lower = unname(exp(ci[1])), upper = unname(exp(ci[2])), p = unname(s[4]), dispersion = sum(residuals(f, "pearson")^2) / f$df.residual) }

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE); d <- if (length(args) >= 1) read.csv(args[1], stringsAsFactors = FALSE) else { cat("Practice data: invented surveillance table\n"); simulate_table() }
  cat("\nResistance in the latest year with exact 95% intervals\n"); ly <- d[d$year == max(d$year), ]
  r <- t(mapply(prop_ci, ly$resistant, ly$tested)); print(data.frame(organism = ly$organism, drug = ly$drug, tested = ly$tested, round(r, 3)), row.names = FALSE)
  cat("\nTrend over years for each organism-drug pair (binomial regression): odds ratio per year\n")
  tr <- do.call(rbind, lapply(split(d, list(d$organism, d$drug), drop = TRUE), function(x) data.frame(pair = paste(x$organism[1], x$drug[1]), t(round(trend_fit(x), 3))))); tr$p_adj <- round(p.adjust(tr$p, "BH"), 4); print(tr, row.names = FALSE)
  cat("\nComparison of two years for one pair (first and last year, Org_1 Drug_C)\n"); x <- d[d$organism == d$organism[1] & d$drug == "Drug_C", ]; x <- x[x$year %in% c(min(x$year), max(x$year)), ]
  pt <- prop.test(x$resistant, x$tested); cat(sprintf("%.3f in %d vs %.3f in %d; difference test p = %.4f\n", x$resistant[1] / x$tested[1], x$year[1], x$resistant[2] / x$tested[2], x$year[2], pt$p.value))
  cat("\nWidth of the interval depends on the number tested: 10% resistance at n = 20, 60, 300, 1000\n"); for (n in c(20, 60, 300, 1000)) { ci <- prop_ci(round(0.1 * n), n); cat(sprintf("n = %4d: %.3f (%.3f to %.3f)\n", n, ci[1], ci[2], ci[3])) }
  cat("\nSurveillance bias: true resistance 10% in all infections, but isolates are mostly cultured after treatment failure\n")
  N <- 1e5; res <- runif(N) < 0.10; p_test <- ifelse(res, 0.60, 0.15); tested <- runif(N) < p_test
  cat(sprintf("Resistance among all infections: %.3f; among the isolates that reached the lab: %.3f (%d isolates)\n", mean(res), mean(res[tested]), sum(tested)))
}
