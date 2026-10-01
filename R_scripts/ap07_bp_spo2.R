# AP-07  Blood pressure and oxygen saturation: repeatability, posture, and agreement between devices. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript ap07_bp_spo2.R [readings.csv]
# readings.csv: subject (anonymous code), device (A or B), posture (seated or standing), rep, sbp, dbp, hr, spo2
# For teaching only. Values are not medical advice and must not be interpreted as anyone's health status.
set.seed(7)
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) d <- read.csv(args[1]) else {
    n <- 16; true_s <- rnorm(n, 112, 9); mk <- function(dev, bias) do.call(rbind, lapply(seq_len(n), function(i) do.call(rbind, lapply(c("seated", "standing"), function(p) data.frame(subject = paste0("S", i), device = dev, posture = p, rep = 1:3,
      sbp = round(true_s[i] + bias + ifelse(p == "standing", 2, 0) + rnorm(3, 0, 4)), dbp = round(true_s[i] * 0.65 + rnorm(3, 0, 3)), hr = round(72 + ifelse(p == "standing", 10, 0) + rnorm(3, 0, 3)), spo2 = pmin(100, round(97.5 + rnorm(3, 0, 0.8))))))))
    d <- rbind(mk("A", 0), mk("B", 4)); cat("Practice data: 16 invented students, 2 devices, seated and standing, 3 readings each. Effects and biases are invented.\n\n") }
  s <- d[d$posture == "seated", ]
  rep_sd <- sqrt(mean(tapply(s$sbp[s$device == "A"], s$subject[s$device == "A"], var))); cat(sprintf("Within-subject repeatability of systolic pressure, device A, seated: sd %.1f mmHg; repeatability coefficient (1.96 x sqrt2 x sd) %.1f mmHg\n", rep_sd, 1.96 * sqrt(2) * rep_sd))
  cat(sprintf("Oxygen saturation, seated, device A: mean %.1f%%, within-subject sd %.2f%%; range of individual means %.1f to %.1f%%\n", mean(s$spo2[s$device == "A"]), sqrt(mean(tapply(s$spo2[s$device == "A"], s$subject[s$device == "A"], var))), min(tapply(s$spo2[s$device == "A"], s$subject[s$device == "A"], mean)), max(tapply(s$spo2[s$device == "A"], s$subject[s$device == "A"], mean))))
  m <- aggregate(cbind(sbp, hr) ~ subject + device + posture, d, mean); a <- m[m$device == "A", ]; w <- reshape(a[, c("subject", "posture", "sbp", "hr")], idvar = "subject", timevar = "posture", direction = "wide")
  for (v in c("sbp", "hr")) { t <- t.test(w[[paste0(v, ".standing")]], w[[paste0(v, ".seated")]], paired = TRUE); cat(sprintf("Standing minus seated %s: %+.1f (95%% CI %.1f to %.1f), p = %.3g\n", v, t$estimate, t$conf.int[1], t$conf.int[2], t$p.value)) }
  ss <- m[m$posture == "seated", ]; dv <- reshape(ss[, c("subject", "device", "sbp")], idvar = "subject", timevar = "device", direction = "wide"); df <- dv$sbp.B - dv$sbp.A
  cat(sprintf("\nAgreement between devices (Bland-Altman, seated systolic): mean difference %.1f mmHg, 95%% limits of agreement %.1f to %.1f; correlation %.2f\n", mean(df), mean(df) - 1.96 * sd(df), mean(df) + 1.96 * sd(df), cor(dv$sbp.A, dv$sbp.B)))
  cat("A high correlation can hide a constant offset: the devices agree in ranking but not in value.\n")
  cat("\nHow many readings? Standard error of a person's mean systolic pressure by number of readings (from the within-subject sd):\n"); for (n in c(1, 2, 3, 5)) cat(sprintf("%d reading(s): %.1f mmHg\n", n, rep_sd / sqrt(n)))
}
