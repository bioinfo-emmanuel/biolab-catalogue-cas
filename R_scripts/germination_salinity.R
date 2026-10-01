# BOT-04  Seed germination and salinity tolerance. Base R only: no packages to install.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
#
# YOUR DATA: a CSV with one row per dish and these columns
#   dish, nacl_mM, seeds, germinated
#   optional daily cumulative counts named d1, d2, d3 ... (germinated seeds counted by that day)
#   nacl_mM     salt concentration in millimolar (0 for the control)
#   seeds       seeds put in the dish
#   germinated  seeds germinated at the end of scoring
# The dish, not the seed, is the replicate.
#
# Run in RStudio: set data_file below and click Source. From a terminal: Rscript germination_salinity.R mydata.csv
# With no file the script makes PRACTICE data (simulated). Never report them as results.

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1) args[1] else NA

if (is.na(data_file)) {
  message("No data file given: using simulated PRACTICE data.")
  set.seed(4)
  conc <- rep(c(0, 25, 50, 100, 150, 200), each = 3)
  p <- plogis(2.2 - 0.022 * conc)                      # true EC50 = 2.2 / 0.022 = 100 mM
  g <- rbinom(length(conc), 10, p)
  daily <- sapply(1:6, function(day) rbinom(length(g), g, pgamma(day, shape = 3, rate = 1.1)))
  daily <- t(apply(daily, 1, cummax)); daily[, 6] <- g
  d <- data.frame(dish = seq_along(conc), nacl_mM = conc, seeds = 10, germinated = g)
  d[, paste0("d", 1:6)] <- daily
} else {
  d <- read.csv(data_file)
}
d$prop <- d$germinated / d$seeds

# ---- 1. logistic (binomial) regression ---------------------------------------------
fit <- glm(cbind(germinated, seeds - germinated) ~ nacl_mM, family = binomial, data = d)
cat("\n1. Binomial regression of germination on salt concentration\n")
print(round(summary(fit)$coefficients, 4))
b0 <- coef(fit)[1]; b1 <- coef(fit)[2]
if (b1 >= 0) stop("Germination did not fall with salt in these data, so EC50 is not defined.")

# ---- 2. EC50: the salt concentration that halves germination (logit(p) = 0) ---------
ec50 <- unname(-b0 / b1)
V <- vcov(fit)
se <- sqrt(unname((V[1, 1] + 2 * ec50 * V[1, 2] + ec50^2 * V[2, 2]) / b1^2))     # delta method
cat(sprintf("\n2. EC50 = %.0f mM (approximate 95%% interval %.0f to %.0f mM)\n", ec50, ec50 - 1.96 * se, ec50 + 1.96 * se))

# ---- 3. is there more scatter than the model allows? -----------------------------------
disp <- deviance(fit) / df.residual(fit)
cat(sprintf("\n3. Residual deviance / df = %.2f. Values well above 1 mean dishes differ more than a simple binomial allows,\n   so treat the interval above as too narrow.\n", disp))

# ---- 4. time to 50% of final germination, if daily counts are present ----------------
dcols <- grep("^d[0-9]+$", names(d), value = TRUE)
if (length(dcols) >= 2) {
  days <- as.numeric(sub("^d", "", dcols))
  t50 <- apply(d[, dcols, drop = FALSE], 1, function(x) {
    if (max(x) == 0) return(NA)
    approx(c(0, x), c(0, days), xout = 0.5 * max(x), ties = "ordered")$y })
  d$t50_days <- t50
  cat("\n4. Time to 50% of final germination (days), mean by concentration\n")
  print(round(tapply(d$t50_days, d$nacl_mM, mean, na.rm = TRUE), 2))
}

# ---- 5. figure ---------------------------------------------------------------------------
x <- seq(0, max(d$nacl_mM), length.out = 100)
plot(d$nacl_mM, d$prop, pch = 19, ylim = c(0, 1), xlab = "NaCl (mM)", ylab = "Proportion germinated")
lines(x, predict(fit, data.frame(nacl_mM = x), type = "response"), col = "red")
abline(v = ec50, lty = 2)
