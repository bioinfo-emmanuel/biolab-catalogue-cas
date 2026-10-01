# BCH-05  Predicting isoelectric points from sequence. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript bch05_isoelectric_point.R [file.fasta]
# Net charge = sum over ionizable groups of the fraction protonated (bases) or deprotonated (acids),
# using the Henderson-Hasselbalch relation. No electrostatic interactions between groups are modelled.

# pKa sets. EMBOSS: Rice, Longden & Bleasby (2000) Trends Genet 16:276-277.
# Grimsley: mean measured values in folded proteins, Grimsley, Gray & Pace (2009) Protein Sci 18:247-251
# (no measured mean for Arg is given there, so Arg uses the EMBOSS value).
pka_sets <- list(
  EMBOSS   = c(Nterm = 8.6, Cterm = 3.6, K = 10.8, R = 12.5, H = 6.5, D = 3.9, E = 4.1, C = 8.5, Y = 10.1),
  Grimsley = c(Nterm = 7.7, Cterm = 3.3, K = 10.5, R = 12.5, H = 6.6, D = 3.5, E = 4.2, C = 6.8, Y = 10.3))

read_fasta <- function(file) { l <- readLines(file); h <- grep("^>", l)
  ends <- c(h[-1] - 1, length(l)); setNames(mapply(function(s, e) paste(l[(s + 1):e], collapse = ""), h, ends), sub("^>", "", l[h])) }

composition <- function(seq) { ch <- strsplit(toupper(gsub("[^A-Za-z]", "", seq)), "")[[1]]
  table(factor(ch, levels = c("K", "R", "H", "D", "E", "C", "Y"))) }

net_charge <- function(pH, comp, pka) {
  pos <- function(pk) 1 / (1 + 10^(pH - pk)); neg <- function(pk) 1 / (1 + 10^(pk - pH))
  pos(pka["Nterm"]) + comp["K"] * pos(pka["K"]) + comp["R"] * pos(pka["R"]) + comp["H"] * pos(pka["H"]) -
    neg(pka["Cterm"]) - comp["D"] * neg(pka["D"]) - comp["E"] * neg(pka["E"]) - comp["C"] * neg(pka["C"]) - comp["Y"] * neg(pka["Y"])
}
pI <- function(seq, set = "EMBOSS") { comp <- composition(seq); pka <- pka_sets[[set]]
  unname(uniroot(function(p) net_charge(p, comp, pka), c(0, 14), tol = 1e-8)$root) }

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  # Check against the published EMBOSS value quoted in MathWorks documentation for this peptide.
  test <- "PQGGGGWGQPHGGGWGQPHGGGGWGQGGSHSQG"
  cat(sprintf("Check peptide: pI = %.4f (documentation value with the EMBOSS pK file: 7.8109)\n\n", pI(test, "EMBOSS")))
  if (length(args) >= 1) { s <- read_fasta(args[1])
    res <- data.frame(name = names(s), length = nchar(s), pI_EMBOSS = sapply(s, pI, "EMBOSS"), pI_Grimsley = sapply(s, pI, "Grimsley"), row.names = NULL)
    print(res, digits = 3) } else {
    # Built-in demo: poly-residue peptides with obvious answers.
    demo <- c(acidic = "DDDDEEEEDDDDEEEE", basic = "KKKKRRRRKKKKRRRR", mixed = "AAKAAEAAHAAYAACAA")
    print(data.frame(peptide = names(demo), pI_EMBOSS = sapply(demo, pI, "EMBOSS"), pI_Grimsley = sapply(demo, pI, "Grimsley"), row.names = NULL), digits = 3) }
  cat("\nNet charge of the check peptide at pH 7.3 (documentation value 0.3629):",
      round(net_charge(7.3, composition(test), pka_sets$EMBOSS), 4), "\n")
}
