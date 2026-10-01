# AP-02  Touch acuity mapping: two-point discrimination across body regions. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript ap02_two_point.R [make_sheet | trials.csv]
#   Rscript ap02_two_point.R make_sheet  prints a randomized presentation order with catch trials.
# trials.csv: subject, region, spacing_mm, stimulus (two or one), response (two or one)
set.seed(2)
make_sheet <- function(regions, spacings, reps = 2) {
  t <- expand.grid(region = regions, spacing_mm = spacings, rep = seq_len(reps), stringsAsFactors = FALSE); t$stimulus <- "two"
  c1 <- data.frame(region = regions, spacing_mm = 0, rep = 1, stimulus = "one", stringsAsFactors = FALSE); c1 <- c1[rep(seq_len(nrow(c1)), each = 3), ]
  s <- rbind(t, c1); s[sample(nrow(s)), ] }
threshold_by_logit <- function(x) {                     # spacing at which two points are reported 50% of the time (two-point trials only)
  x <- x[x$stimulus == "two", ]; x$y <- as.integer(x$response == "two"); if (length(unique(x$y)) < 2) return(NA)
  f <- glm(y ~ spacing_mm, binomial, x); b <- coef(f); if (b[2] <= 0) NA else -b[1] / b[2] }
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1 && args[1] == "make_sheet") { s <- make_sheet(c("fingertip", "palm", "forearm", "back_of_hand"), c(2, 4, 8, 15, 30, 50)); print(head(s, 12), row.names = FALSE); cat("... total presentations:", nrow(s), "(including catch trials with one point)\n"); quit(save = "no") }
  if (length(args) >= 1) d <- read.csv(args[1]) else {
    reg <- c(fingertip = 3, palm = 10, back_of_hand = 22, forearm = 35); n <- 10; d <- do.call(rbind, lapply(seq_len(n), function(i) { s <- make_sheet(names(reg), c(2, 4, 8, 15, 30, 50), 2); s$subject <- paste0("S", i)
      s$response <- mapply(function(r, sp, st) if (st == "one") ifelse(runif(1) < 0.08, "two", "one") else ifelse(runif(1) < plogis((sp - reg[[r]] * rnorm(1, 1, 0.15)) / (0.25 * reg[[r]]) ), "two", "one"), s$region, s$spacing_mm, s$stimulus); s }))
    cat("Practice data: 10 invented subjects, 4 regions, 6 spacings x 2 repeats plus catch trials each; region thresholds are invented\n\n") }
  cat(sprintf("False-alarm rate on one-point catch trials: %.1f%% (a high value means the subject reports two points too readily)\n\n", 100 * mean(d$response[d$stimulus == "one"] == "two")))
  th <- do.call(rbind, lapply(split(d, list(d$subject, d$region)), function(x) data.frame(subject = x$subject[1], region = x$region[1], threshold_mm = threshold_by_logit(x)))); th <- th[!is.na(th$threshold_mm), ]
  res <- aggregate(threshold_mm ~ region, th, function(x) c(n = length(x), median = median(x), q1 = quantile(x, .25), q3 = quantile(x, .75))); print(do.call(data.frame, res), digits = 3)
  w <- reshape(th, idvar = "subject", timevar = "region", direction = "wide"); w <- w[complete.cases(w), ]; k <- friedman.test(as.matrix(w[, -1])); cat(sprintf("\nFriedman test of threshold across regions (subjects as blocks): chi-square = %.1f, df = %d, p = %.2e\n", k$statistic, k$parameter, k$p.value))
  cat("Fingertip vs forearm, paired Wilcoxon: p =", signif(wilcox.test(w$threshold_mm.fingertip, w$threshold_mm.forearm, paired = TRUE, exact = FALSE)$p.value, 2), "\n")
  cat("\nBlinding and randomization matter: with the presentation order fixed (all wide spacings first), subjects learn the pattern and the thresholds shift. Use the randomized sheet from make_sheet.\n")
}
