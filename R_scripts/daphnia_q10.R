# ZOO-04  Heart rate and temperature in Daphnia: calculating Q10. Base R only: no packages to install.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
#
# YOUR DATA: a CSV with one row per animal and temperature, and these columns
#   animal, temp_c, beats_1, beats_2, seconds
#   temp_c       water temperature in degrees C (use the average of the readings before and after counting)
#   beats_1, beats_2   two counts of heartbeats over `seconds`
#   seconds      length of each count (for example 15)
# The animal is the block: every animal is measured at several temperatures in a random order.
#
# Run in RStudio: set data_file below and click Source. From a terminal: Rscript daphnia_q10.R mydata.csv
# With no file the script makes PRACTICE data (simulated). Never report them as results.

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1) args[1] else NA

if (is.na(data_file)) {
  message("No data file given: using simulated PRACTICE data.")
  set.seed(9)
  temps <- c(10, 15, 20, 25, 30); animals <- 8; q10 <- 2.2
  d <- do.call(rbind, lapply(1:animals, function(a) { base <- rnorm(1, 200, 18)                                  # beats per minute at 20 C
    data.frame(animal = a, temp_c = temps + rnorm(5, 0, 0.3), seconds = 15,
               beats_1 = rpois(5, base * q10^((temps - 20) / 10) / 4), beats_2 = rpois(5, base * q10^((temps - 20) / 10) / 4)) }))
} else {
  d <- read.csv(data_file)
}
d$rate <- (d$beats_1 + d$beats_2) / 2 / d$seconds * 60                   # beats per minute
d$animal <- factor(d$animal)

cat("\n1. Mean heart rate (beats per minute) by nominal temperature\n")
d$T_group <- round(d$temp_c / 5) * 5
print(round(tapply(d$rate, d$T_group, mean), 1))

cat("\n2. Q10 from all the data: regress ln(rate) on temperature with animal as a block\n")
fit <- lm(log(rate) ~ temp_c + animal, data = d)
b <- coef(fit)[["temp_c"]]; ci <- confint(fit)["temp_c", ]
cat(sprintf("slope = %.4f per degree C; Q10 = exp(10 x slope) = %.2f (95%% interval %.2f to %.2f)\n", b, exp(10 * b), exp(10 * ci[1]), exp(10 * ci[2])))

cat("\n3. Q10 between neighbouring temperatures (mean rate at each): a sudden fall flags stress\n")
mr <- tapply(d$rate, d$T_group, mean); tg <- as.numeric(names(mr))
q <- (mr[-1] / mr[-length(mr)])^(10 / diff(tg)); names(q) <- paste0(tg[-length(tg)], "-", tg[-1], " C"); print(round(q, 2))

cat("\n4. Q10 from two temperatures, worked example: 150 bpm at 15 C and 300 bpm at 25 C ->", round((300 / 150)^(10 / (25 - 15)), 2), "\n")
plot(d$temp_c, log(d$rate), pch = 19, xlab = "Temperature (C)", ylab = "ln(heart rate, beats per minute)")
abline(lm(log(rate) ~ temp_c, data = d), col = "red")
