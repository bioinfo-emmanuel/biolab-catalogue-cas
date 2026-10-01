# ACH-02  Acetic acid in vinegar by titration, with standardization of NaOH. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript ach02_vinegar_titration.R [standardization.csv vinegar.csv]
# standardization.csv: khp_g, naoh_ml    (mass of primary standard potassium hydrogen phthalate, volume of NaOH to endpoint)
# vinegar.csv: brand, label_pct, aliquot_ml, dilution_factor, naoh_ml   (label_pct = % acetic acid stated on the bottle)
set.seed(202)
M_KHP <- 204.22; M_HAC <- 60.05          # g/mol (KHP is C8H5KO4; acetic acid is C2H4O2)

standardize <- function(s) { m <- (s$khp_g / M_KHP) / (s$naoh_ml / 1000); c(mean = mean(m), sd = sd(m), n = length(m), rsd_pct = 100 * sd(m) / mean(m)) }
# g of acetic acid per 100 mL of undiluted vinegar (about equal to % w/v; % w/w needs the density)
acid_g_per_100ml <- function(naoh_ml, M, aliquot_ml, dil) M * (naoh_ml / 1000) * M_HAC / (aliquot_ml / dil) * 100 / 1

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 2) { st <- read.csv(args[1]); vin <- read.csv(args[2]) } else {
    trueM <- 0.1012; khp <- round(runif(4, 0.40, 0.55), 4); st <- data.frame(khp_g = khp, naoh_ml = round(khp / M_KHP / trueM * 1000 + rnorm(4, 0, 0.05), 2))
    real <- c(BrandA = 4.9, BrandB = 5.6, BrandC = 3.6); lab <- c(BrandA = 5, BrandB = 5, BrandC = 4)
    # NaOH volume (mL) for a 10 mL aliquot of the 1-in-10 dilution: acid in aliquot = real/100 g/mL x (10/10) mL
    vol <- function(pct) pct / 100 * (10 / 10) / M_HAC / trueM * 1000
    vin <- do.call(rbind, lapply(names(real), function(b) data.frame(brand = b, label_pct = lab[[b]], aliquot_ml = 10, dilution_factor = 10, naoh_ml = round(vol(real[[b]]) + rnorm(4, 0, 0.08), 2))))
    cat("Practice data: 4 standardization runs and 4 titrations for each of 3 vinegars (a 10 mL sample diluted to 100 mL, 10 mL aliquots)\n\n") }
  s <- standardize(st); cat(sprintf("NaOH standardization: %.4f M (sd %.4f, rsd %.2f%%, n = %d)\n\n", s["mean"], s["sd"], s["rsd_pct"], s["n"]))
  vin$acid <- acid_g_per_100ml(vin$naoh_ml, s["mean"], vin$aliquot_ml, vin$dilution_factor)
  res <- do.call(rbind, lapply(split(vin, vin$brand), function(x) { t <- t.test(x$acid, mu = x$label_pct[1]); data.frame(brand = x$brand[1], label = x$label_pct[1], n = nrow(x), mean = mean(x$acid), sd = sd(x$acid), ci_low = t$conf.int[1], ci_high = t$conf.int[2], p_vs_label = t$p.value) }))
  print(format(res, digits = 3), row.names = FALSE)
  cat("\nWhat if the NaOH is not standardized and its concentration is assumed to be exactly 0.1000 M (true value about", round(s["mean"], 4), "M)?\n")
  print(round(tapply(acid_g_per_100ml(vin$naoh_ml, 0.1, vin$aliquot_ml, vin$dilution_factor), vin$brand, mean), 2))
}
