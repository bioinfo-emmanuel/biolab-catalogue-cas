# ZOO-06 and CVA-05  Metabolic and trait scaling: log-log regression, exponent tests, and outliers. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript zoo06_allometry.R [data.csv]
# data.csv columns: species, mass, trait   (positive numbers; optionally a column group).
# The same script serves ZOO-06 (trait = metabolic rate) and CVA-05 (trait = brain, heart or bone dimension).
# With no argument it simulates 60 species with a true exponent of 0.75 and lognormal scatter.
set.seed(75)

fit_scaling <- function(d) {
  d$lm <- log10(d$mass); d$lt <- log10(d$trait)
  ols <- lm(lt ~ lm, d); r <- cor(d$lm, d$lt)
  sma <- sign(r) * sd(d$lt) / sd(d$lm)                          # standardized major axis slope
  ci <- confint(ols)["lm", ]
  list(ols = ols, ols_slope = coef(ols)[["lm"]], ci = ci, sma_slope = sma, r2 = r^2,
       # SMA slope confidence interval (Pitman 1939): uses r and n
       sma_ci = { n <- nrow(d); B <- (qt(0.975, n - 2)^2 * (1 - r^2)) / (n - 2); sma * (sqrt(B + 1) + c(-1, 1) * sqrt(B)) })
}
test_exponent <- function(fit, b0) { se <- summary(fit$ols)$coefficients["lm", 2]; df <- fit$ols$df.residual
  t <- (fit$ols_slope - b0) / se; c(hypothesis = b0, t = t, p = 2 * pt(-abs(t), df)) }
outliers <- function(d, fit, k = 2) { s <- rstandard(fit$ols); d[abs(s) > k, ] }

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) d <- read.csv(args[1]) else {
    n <- 60; mass <- 10^runif(n, -1, 4)
    d <- data.frame(species = paste0("sp", 1:n), mass = mass, trait = 3.5 * mass^0.75 * 10^rnorm(n, 0, 0.12))
    d$trait[7] <- d$trait[7] * 3; cat("Practice data: 60 simulated species, true exponent 0.75, one species inflated threefold\n\n") }
  f <- fit_scaling(d)
  cat(sprintf("n = %d  R squared = %.3f\n", nrow(d), f$r2))
  cat(sprintf("OLS slope = %.3f  (95%% CI %.3f to %.3f)\nSMA slope = %.3f  (95%% CI %.3f to %.3f)\n", f$ols_slope, f$ci[1], f$ci[2], f$sma_slope, f$sma_ci[1], f$sma_ci[2]))
  cat("\nTests of the slope against two expected exponents (OLS):\n"); print(round(rbind(test_exponent(f, 0.75), test_exponent(f, 2/3)), 4))
  cat("\nStandardized residuals beyond 2:\n"); print(outliers(d, f))
  cat("\nMeasurement error in mass (log10 units) pulls the OLS slope toward zero (regression dilution); the effect only shows when the error is large relative to the spread of masses\n")
  for (err in c(0, 0.5, 1)) { m2 <- d$mass * 10^rnorm(nrow(d), 0, err); ff <- fit_scaling(data.frame(mass = m2, trait = d$trait)); cat(sprintf("mass error sd %.1f (log10 units): OLS %.3f  SMA %.3f\n", err, ff$ols_slope, ff$sma_slope)) }
  cat("\nEffect of dropping the outlier: "); f2 <- fit_scaling(d[-7, ]); cat(sprintf("OLS slope %.3f (CI %.3f to %.3f)\n", f2$ols_slope, f2$ci[1], f2$ci[2]))
  png("scaling_plot.png", 700, 500); plot(log10(d$mass), log10(d$trait), xlab = "log10 mass", ylab = "log10 trait"); abline(f$ols, col = "red"); dev.off()
}
