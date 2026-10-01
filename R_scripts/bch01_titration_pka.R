# BCH-01  Titration of a weak acid: equivalence point, pKa, and buffer capacity. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript bch01_titration_pka.R [titration.csv]
# titration.csv: naoh_ml, ph   (volume of NaOH added and the pH after each addition)
# Optional: set acid_ml, base_M and the acid concentration below for the buffer-capacity calculation.
set.seed(11)
acid_ml <- 25; base_M <- 0.100

# Model titration of a monoprotic weak acid (charge balance solved for [H+]); used to make practice data and to check the analysis.
model_ph <- function(vb_ml, ca_M, va_ml, cb_M, pKa, Kw = 1e-14) {
  Ka <- 10^-pKa; vt <- va_ml + vb_ml; Ca <- ca_M * va_ml / vt; Cb <- cb_M * vb_ml / vt
  sapply(seq_along(vt), function(i) { f <- function(lh) { h <- 10^lh; Cb[i] + h - Kw / h - Ca[i] * Ka / (Ka + h) }; uniroot(f, c(-13, 0))$root * -1 }) }

analyze <- function(d) {
  d <- d[order(d$naoh_ml), ]; v <- d$naoh_ml; p <- d$ph; n <- length(v)
  slope <- diff(p) / diff(v); vmid <- (v[-1] + v[-n]) / 2
  i <- which.max(slope); ve1 <- vmid[i]                                             # first-derivative estimate of the equivalence volume
  # refine with the second derivative zero crossing between neighbouring midpoints
  s2 <- diff(slope) / diff(vmid); v2 <- (vmid[-1] + vmid[-length(vmid)]) / 2
  j <- which(s2[-length(s2)] > 0 & s2[-1] <= 0 & v2[-length(v2)] > ve1 - 1.5 & v2[-length(v2)] < ve1 + 1.5)[1]
  ve2 <- if (!is.na(j)) approx(s2[j:(j + 1)], v2[j:(j + 1)], xout = 0)$y else NA
  ve <- if (!is.na(ve2)) ve2 else ve1; half <- ve / 2; pKa <- approx(v, p, xout = half)$y
  list(ve_first = ve1, ve_second = ve2, ve = ve, half = half, pKa = pKa, slope = slope, vmid = vmid) }
buffer_capacity <- function(d, va_ml, base_M) { d <- d[order(d$naoh_ml), ]; n <- nrow(d); dn <- diff(d$naoh_ml) / 1000 * base_M; vt <- (va_ml + (d$naoh_ml[-1] + d$naoh_ml[-n]) / 2) / 1000
  data.frame(ph = (d$ph[-1] + d$ph[-n]) / 2, beta = (dn / vt) / diff(d$ph)) }

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) d <- read.csv(args[1]) else {
    ca <- 0.08; va <- acid_ml; ve_true <- ca * va / base_M
    vv <- sort(unique(c(seq(0, 16, by = 1), seq(16, 24, by = 0.5), seq(24, 26, by = 0.2), seq(26, 32, by = 1))))
    d <- data.frame(naoh_ml = vv, ph = round(model_ph(vv, ca, va, base_M, 4.76) + rnorm(length(vv), 0, 0.02), 2))
    cat(sprintf("Practice data: 25 mL of 0.0800 M weak acid (true pKa 4.76) titrated with 0.1000 M NaOH; true equivalence volume %.1f mL; pH read to 0.01 with noise of 0.02\n\n", ve_true)) }
  r <- analyze(d)
  cat(sprintf("Equivalence volume from the steepest point: %.2f mL (second-derivative estimate %s mL)\n", r$ve_first, ifelse(is.na(r$ve_second), "not found", sprintf("%.2f", r$ve_second))))
  cat(sprintf("Half of the equivalence volume (%.2f mL): pH = %.2f, so pKa is about %.2f\n", r$half, r$pKa, r$pKa))
  cat(sprintf("Concentration of the acid from the equivalence volume: %.4f M\n", base_M * r$ve / acid_ml))
  bc <- buffer_capacity(d, acid_ml, base_M); k <- which.max(bc$beta[bc$ph < 8]); cat(sprintf("Buffer capacity peaks at pH %.2f (near the pKa) at about %.3f mol per L per pH unit\n", bc$ph[k], bc$beta[k]))
  cat("\nEffect of a pH meter that is off by 0.30 pH units (no calibration): pKa read as", round(r$pKa + 0.3, 2), "but the equivalence volume is unchanged:", round(r$ve, 2), "mL\n")
  cat("Coarse addition (1 mL steps everywhere) equivalence volume:", round(analyze(d[d$naoh_ml %in% seq(0, 32, 1), ])$ve, 1), "mL\n")
}
