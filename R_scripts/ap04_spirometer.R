# AP-04  Lung volumes with a homemade water-displacement spirometer. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript ap04_spirometer.R [calibration.csv volumes.csv]
# calibration.csv: added_ml, mark_cm   (volume of water poured in with a cylinder, and the water-level mark on the bottle)
# volumes.csv: subject (anonymous code), sex, height_cm, measure (tidal, erv or vc), trial, ml
set.seed(4)
calibrate <- function(cal) { f <- lm(added_ml ~ mark_cm, cal); list(fit = f, to_ml = function(mark_cm) predict(f, data.frame(mark_cm = mark_cm))) }
best_two_agree <- function(x, tol = 0.05) { x <- sort(x, decreasing = TRUE); if (length(x) >= 2 && (x[1] - x[2]) / x[1] <= tol) mean(x[1:2]) else NA }

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 2) { cal <- read.csv(args[1]); v <- read.csv(args[2]) } else {
    cal <- data.frame(added_ml = seq(0, 4500, by = 500), mark_cm = round(seq(0, 30, length.out = 10) + rnorm(10, 0, 0.1), 1)); cal$added_ml <- seq(0, 4500, by = 500)
    n <- 24; sex <- rep(c("F", "M"), each = n / 2); h <- round(ifelse(sex == "F", rnorm(n, 157, 6), rnorm(n, 169, 6)))
    vc <- round(ifelse(sex == "F", -1700 + 32 * h, -2400 + 39 * h) * rnorm(n, 1, 0.07) / 10) * 10
    v <- do.call(rbind, lapply(seq_len(n), function(i) rbind(data.frame(subject = paste0("S", i), sex = sex[i], height_cm = h[i], measure = "vc", trial = 1:3, ml = round(vc[i] * c(1, 0.99, 0.96) * rnorm(3, 1, 0.015) / 50) * 50),
      data.frame(subject = paste0("S", i), sex = sex[i], height_cm = h[i], measure = "tidal", trial = 1:3, ml = round(rnorm(3, 500, 60) / 50) * 50), data.frame(subject = paste0("S", i), sex = sex[i], height_cm = h[i], measure = "erv", trial = 1:3, ml = round(rnorm(3, ifelse(sex[i] == "F", 900, 1200), 150) / 50) * 50))))
    cat("Practice data: 24 invented students (12 F, 12 M), three trials of each measure, exhaled volumes to the nearest 50 mL. The relationships are invented for practice, not reference values.\n\n") }
  if ("mark_cm" %in% names(cal)) { k <- calibrate(cal); cat(sprintf("Calibration: volume = %.1f + %.1f x mark (cm); R squared %.5f; residual sd %.1f mL\n", coef(k$fit)[1], coef(k$fit)[2], summary(k$fit)$r.squared, sigma(k$fit))) }
  agg <- aggregate(ml ~ subject + sex + height_cm + measure, v, function(x) { y <- best_two_agree(x); if (is.na(y)) mean(x) else y }); cat(sprintf("Repeatability: share of subjects whose two largest vital-capacity trials agree within 5%%: %.0f%%\n\n", 100 * mean(tapply(v$ml[v$measure == "vc"], v$subject[v$measure == "vc"], function(x) !is.na(best_two_agree(x))))))
  w <- reshape(agg, idvar = c("subject", "sex", "height_cm"), timevar = "measure", direction = "wide", v.names = "ml")
  ag <- do.call(data.frame, aggregate(cbind(ml.vc, ml.tidal, ml.erv) ~ sex, w, function(x) c(mean = mean(x), sd = sd(x)))); ag[, -1] <- round(ag[, -1], 0); print(ag, row.names = FALSE)
  f <- lm(ml.vc ~ height_cm + sex, w); cat("\nVital capacity regressed on height and sex:\n"); print(round(summary(f)$coefficients, 3)); cat(sprintf("R squared %.2f; residual sd %.0f mL\n", summary(f)$r.squared, sigma(f)))
  cat(sprintf("Each extra cm of height adds about %.0f mL of vital capacity (95%% CI %.0f to %.0f); males are %.0f mL higher at the same height\n", coef(f)["height_cm"], confint(f)["height_cm", 1], confint(f)["height_cm", 2], coef(f)["sexM"]))
  cat("Bottle reading resolution: a 1 mm mark error on a 30 cm scale for 4500 mL is about", round(4500 / 300, 0), "mL per mm.\n")
}
