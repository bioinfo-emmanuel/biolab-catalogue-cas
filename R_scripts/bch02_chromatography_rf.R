# BCH-02  Amino acids by paper chromatography: Rf values and matching unknowns. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript bch02_chromatography_rf.R [spots.csv]
# spots.csv: run, sample, spot_cm, front_cm   (distance from the origin line to the centre of the spot, and to the solvent front)
# Standards are named by amino acid; unknown samples start with "unk". With no argument the script simulates a run.
set.seed(22)
rf <- function(spot, front) spot / front
match_unknown <- function(rf_unknown, std, tol = 0.05) { d <- abs(std$rf_mean - rf_unknown); ok <- which(d <= tol); if (length(ok)) std$sample[ok[order(d[ok])]] else NA }

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) d <- read.csv(args[1]) else {
    true_rf <- c(leucine = 0.72, alanine = 0.30, glycine = 0.21, serine = 0.17, aspartic_acid = 0.10)
    mk <- function(run, s, r) { fr <- runif(1, 11.5, 12.5); data.frame(run = run, sample = s, spot_cm = round(r * fr + rnorm(1, 0, 0.12), 2), front_cm = round(fr, 2)) }
    d <- do.call(rbind, c(lapply(1:4, function(k) do.call(rbind, lapply(names(true_rf), function(s) mk(k, s, true_rf[[s]])))),
                          lapply(1:4, function(k) mk(k, "unk_1", true_rf[["alanine"]])), lapply(1:4, function(k) mk(k, "unk_2", true_rf[["leucine"]]))))
    cat("Practice data: 4 runs of 5 standards and 2 unknowns. The Rf values are invented for practice; measure your own standards on your own paper and solvent.\n\n") }
  d$rf <- rf(d$spot_cm, d$front_cm)
  st <- d[!grepl("^unk", d$sample), ]; un <- d[grepl("^unk", d$sample), ]
  tab <- do.call(rbind, lapply(split(st, st$sample), function(x) data.frame(sample = x$sample[1], n = nrow(x), rf_mean = mean(x$rf), rf_sd = sd(x$rf))))
  tab <- tab[order(-tab$rf_mean), ]; print(format(tab, digits = 3), row.names = FALSE)
  cat("\nUnknowns (mean Rf, and the standard(s) within +/- 0.05):\n")
  for (u in unique(un$sample)) { m <- mean(un$rf[un$sample == u]); cat(sprintf("%s: Rf %.3f -> %s\n", u, m, paste(match_unknown(m, tab), collapse = ", "))) }
  cat("\nWithin-run comparison (each run's standard leucine against that run's unknown_2): Rf differences by run\n")
  for (k in unique(d$run)) { a <- d$rf[d$run == k & d$sample == "leucine"]; b <- d$rf[d$run == k & d$sample == "unk_2"]; if (length(a) && length(b)) cat(sprintf("run %d: %.3f\n", k, b - a)) }
  cat("\nRf from different runs varies with solvent front position, temperature and paper; always run standards on the same sheet as the unknown.\n")
  cat(sprintf("Largest between-run sd of a standard: %.3f\n", max(tab$rf_sd)))
}
