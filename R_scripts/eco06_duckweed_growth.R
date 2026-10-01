# ECO-06  Population growth in duckweed: exponential and logistic models. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript eco06_duckweed_growth.R [counts.csv]
# counts.csv: treatment, container, day, fronds
set.seed(6)
logistic <- function(t, N0, K, r) K / (1 + ((K - N0) / N0) * exp(-r * t))
fit_container <- function(x) {
  x <- x[order(x$day), ]; t <- x$day - min(x$day); N <- x$fronds
  early <- N <= 0.5 * max(N); e <- if (sum(early) >= 3) lm(log(N[early]) ~ t[early]) else lm(log(N) ~ t)
  r0 <- max(coef(e)[2], 0.02); K0 <- max(N) * 1.1
  lg <- tryCatch(nls(N ~ logistic(t, N0, K, r), start = list(N0 = N[1], K = K0, r = r0), control = nls.control(maxiter = 200, warnOnly = TRUE)), error = function(e) NULL)
  exp_fit <- lm(log(N) ~ t); n <- length(N); aic_exp <- n * log(sum((N - exp(fitted(exp_fit)))^2) / n) + 2 * 2
  c(r_early = unname(coef(e)[2]), doubling_days = log(2) / unname(coef(e)[2]), K = if (!is.null(lg)) coef(lg)[["K"]] else NA, r_logistic = if (!is.null(lg)) coef(lg)[["r"]] else NA, aic_logistic = if (!is.null(lg)) AIC(lg) else NA, aic_exponential = AIC(nls(N ~ a * exp(b * t), start = list(a = N[1], b = coef(exp_fit)[[2]])))) }

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) d <- read.csv(args[1]) else {
    days <- seq(0, 27, by = 3); pars <- list(plain = c(r = 0.20, K = 180), fertilised = c(r = 0.26, K = 420))
    d <- do.call(rbind, lapply(names(pars), function(tr) do.call(rbind, lapply(1:4, function(cn) { p <- pars[[tr]]; N <- logistic(days, 6, p[["K"]] * rnorm(1, 1, 0.08), p[["r"]] * rnorm(1, 1, 0.06)); data.frame(treatment = tr, container = paste0(tr, cn), day = days, fronds = pmax(1, round(N * rnorm(length(days), 1, 0.05)))) }))))
    cat("Practice data: 4 containers in each of 2 invented treatments, starting with 6 fronds and counted every 3 days for 27 days\n\n") }
  res <- do.call(rbind, lapply(split(d, d$container), function(x) { r <- fit_container(x); data.frame(treatment = x$treatment[1], container = x$container[1], t(round(r, 3))) })); print(res, row.names = FALSE)
  cat("\nTreatment means:\n"); tm <- aggregate(cbind(r_early, doubling_days, K, r_logistic) ~ treatment, res, mean); tm[, -1] <- round(tm[, -1], 3); print(tm, row.names = FALSE)
  for (v in c("r_early", "K")) { t <- t.test(res[[v]] ~ res$treatment); cat(sprintf("Welch t test of %s between treatments: p = %.4f\n", v, t$p.value)) }
  cat(sprintf("\nModel comparison over containers: logistic fits better by AIC in %d of %d containers\n", sum(res$aic_logistic < res$aic_exponential, na.rm = TRUE), nrow(res)))
  cat("Check: N0 = 6, K = 180 and r = 0.2 per day gives", round(logistic(9, 6, 180, 0.2)), "fronds on day 9; the logistic curve reaches K/2 at t = ln((K - N0) / N0) / r =", round(log((180 - 6) / 6) / 0.2, 1), "days\n")
}
