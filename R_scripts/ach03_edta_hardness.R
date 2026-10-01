# ACH-03  Total water hardness by EDTA titration, with spike recovery. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript ach03_edta_hardness.R [titrations.csv]
# titrations.csv: sample, edta_ml, sample_ml, spike_mg_l (0 for unspiked), and a first line comment is not allowed.
# Set the EDTA molarity below (or pass it as a second argument) after standardizing against a calcium standard.
set.seed(303)
M_CACO3 <- 100.09                      # g/mol; hardness is reported as mg/L CaCO3
hardness <- function(edta_ml, sample_ml, M_edta) M_edta * edta_ml * M_CACO3 * 1000 / sample_ml     # mg/L as CaCO3
standardize_edta <- function(caco3_mg, edta_ml) (caco3_mg / M_CACO3) / edta_ml                      # mol/L (mg / (g/mol) = mmol; mmol / mL = mol/L)
recovery <- function(found_spiked, found_unspiked, added) 100 * (found_spiked - found_unspiked) / added

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE); M <- 0.0100
  if (length(args) >= 1) d <- read.csv(args[1]) else {
    truth <- c(tap = 85, well = 210, bottled = 28, stream = 120)
    d <- do.call(rbind, lapply(names(truth), function(s) data.frame(sample = s, sample_ml = 50, spike_mg_l = 0, edta_ml = round(truth[[s]] * 50 / (M * M_CACO3 * 1000) + rnorm(3, 0, 0.06), 2))))
    sp <- data.frame(sample = "tap", sample_ml = 50, spike_mg_l = 40, edta_ml = round((85 + 40 * 0.96) * 50 / (M * M_CACO3 * 1000) + rnorm(3, 0, 0.06), 2)); d <- rbind(d, sp)
    cat("Practice data: 0.0100 M EDTA, 50 mL samples, 3 titrations each; the tap sample is also spiked with 40 mg/L as CaCO3\n\n") }
  d$hardness <- hardness(d$edta_ml, d$sample_ml, M)
  un <- d[d$spike_mg_l == 0, ]
  res <- do.call(rbind, lapply(split(un, un$sample), function(x) { t <- t.test(x$hardness); data.frame(sample = x$sample[1], n = nrow(x), mean = mean(x$hardness), sd = sd(x$hardness), ci_low = t$conf.int[1], ci_high = t$conf.int[2]) }))
  res$class <- cut(res$mean, c(0, 60, 120, 180, Inf), c("soft", "moderately hard", "hard", "very hard"), right = FALSE); print(format(res, digits = 3), row.names = FALSE)
  sp <- d[d$spike_mg_l > 0, ]; if (nrow(sp)) for (s in unique(sp$sample)) { r <- recovery(mean(sp$hardness[sp$sample == s]), mean(un$hardness[un$sample == s]), sp$spike_mg_l[sp$sample == s][1]); cat(sprintf("\nSpike recovery for %s: %.1f%% (acceptable range is commonly taken as roughly 80 to 120%%)\n", s, r)) }
  cat("\nStandardization check: 10.0 mg CaCO3 titrated with 10.0 mL EDTA gives", round(standardize_edta(10, 10), 5), "M\n")
  cat("Effect of a 0.10 mL endpoint error on a 50 mL sample with 0.0100 M EDTA:", round(hardness(0.10, 50, 0.01), 1), "mg/L\n")
  cat("Hardness classes used here (mg/L CaCO3): soft under 60, moderately hard 60 to 120, hard 120 to 180, very hard 180 and above\n")
}
