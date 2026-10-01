# AP-01  Reaction time and sensory processing: the ruler-drop test. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript ap01_reaction_time.R [trials.csv]
# trials.csv: subject (anonymous code), trial, hand (dominant or non_dominant), condition (baseline or distraction), catch_cm (distance the ruler fell before it was caught)
# Time from distance: t = sqrt(2 d / g), g = 9.81 m/s2, d in metres.
set.seed(1)
g <- 9.81
rt_s <- function(catch_cm) sqrt(2 * catch_cm / 100 / g)
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) d <- read.csv(args[1]) else {
    n <- 12; subj <- paste0("S", 1:n); base <- rnorm(n, 0.240, 0.020)
    mk <- function(s, i, hand, cond) { ti <- 1:20; t <- base[i] + ifelse(hand == "non_dominant", 0.012, 0) + ifelse(cond == "distraction", 0.030, 0) - 0.06 * (1 - exp(-ti / 6)) + rnorm(20, 0, 0.018); data.frame(subject = s, trial = ti, hand = hand, condition = cond, catch_cm = round(pmax(t, 0.08)^2 * g / 2 * 100, 0)) }
    d <- do.call(rbind, unlist(lapply(seq_len(n), function(i) list(mk(subj[i], i, "dominant", "baseline"), mk(subj[i], i, "non_dominant", "baseline"), mk(subj[i], i, "dominant", "distraction"))), recursive = FALSE))
    cat("Practice data: 12 invented students, 20 trials in each of three conditions, catch distance to the nearest cm. The effects are invented.\n\n") }
  d$rt <- rt_s(d$catch_cm)
  cat("Check: a ruler caught at 15 cm gives", round(rt_s(15) * 1000), "ms; at 20 cm", round(rt_s(20) * 1000), "ms; a 1 cm difference near 20 cm is worth", round((rt_s(21) - rt_s(20)) * 1000, 1), "ms\n\n")
  m <- aggregate(rt ~ subject + hand + condition, d, mean)
  base <- m[m$condition == "baseline", ]; b <- reshape(base[, c("subject", "hand", "rt")], idvar = "subject", timevar = "hand", direction = "wide")
  cat(sprintf("Baseline mean reaction time: dominant %.0f ms, non-dominant %.0f ms\n", 1000 * mean(b$rt.dominant), 1000 * mean(b$rt.non_dominant)))
  tt <- t.test(b$rt.dominant, b$rt.non_dominant, paired = TRUE); cat(sprintf("Paired t test (subject means): difference %.1f ms (95%% CI %.1f to %.1f), p = %.4f\n", 1000 * tt$estimate, 1000 * tt$conf.int[1], 1000 * tt$conf.int[2], tt$p.value))
  dd <- m[m$hand == "dominant", ]; dw <- reshape(dd[, c("subject", "condition", "rt")], idvar = "subject", timevar = "condition", direction = "wide")
  tt2 <- t.test(dw$rt.distraction, dw$rt.baseline, paired = TRUE); cat(sprintf("Distraction vs baseline (dominant hand): +%.1f ms (95%% CI %.1f to %.1f), p = %.2e\n", 1000 * tt2$estimate, 1000 * tt2$conf.int[1], 1000 * tt2$conf.int[2], tt2$p.value))
  x <- d[d$hand == "dominant" & d$condition == "baseline", ]; pm <- tapply(x$rt, x$trial, mean); cat("\nPractice effect (dominant hand, baseline): mean reaction time in ms by trial block\n"); print(round(1000 * tapply(x$rt, cut(x$trial, c(0, 5, 10, 15, 20)), mean), 0))
  f <- lm(rt ~ trial, x); cat(sprintf("Linear trend %.2f ms per trial (p = %.3g)\n", 1000 * coef(f)[2], summary(f)$coefficients[2, 4]))
  cat(sprintf("\nWithin-person variability (sd of trials, dominant baseline): median %.0f ms; between-person sd of means: %.0f ms\n", 1000 * median(tapply(x$rt, x$subject, sd)), 1000 * sd(tapply(x$rt, x$subject, mean))))
  cat("The paired design removes between-person differences; pooling all trials as independent would treat 20 trials from one person as 20 people.\n")
}
