# BIO-04  Designing a cloning experiment in silico: restriction sites, primers, virtual gel. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript bio04_cloning_design.R [insert.fasta plasmid.fasta]
# The plasmid is treated as circular, the insert as linear. With no arguments the script builds SYNTHETIC
# sequences (a 900 bp "gene" and a 3200 bp "plasmid" with a multiple cloning site) so it runs offline.
set.seed(404)

# Recognition sequence and cut position (bases after which the top strand is cut) for common enzymes.
enzymes <- data.frame(name = c("EcoRI", "BamHI", "HindIII", "XhoI", "NdeI", "NcoI", "SalI", "KpnI", "SacI", "XbaI", "NotI", "PstI", "SmaI", "EcoRV"),
  site = c("GAATTC", "GGATCC", "AAGCTT", "CTCGAG", "CATATG", "CCATGG", "GTCGAC", "GGTACC", "GAGCTC", "TCTAGA", "GCGGCCGC", "CTGCAG", "CCCGGG", "GATATC"),
  cut = c(1, 1, 1, 1, 2, 1, 1, 5, 5, 1, 2, 5, 3, 3), stringsAsFactors = FALSE)

read_fasta1 <- function(f) { l <- readLines(f); toupper(paste(l[!grepl("^>", l)], collapse = "")) }
find_sites <- function(seq, site, circular = FALSE) { s <- if (circular) paste0(seq, substr(seq, 1, nchar(site) - 1)) else seq
  m <- gregexpr(site, s, fixed = TRUE)[[1]]; if (m[1] == -1) integer(0) else as.integer(m)[as.integer(m) <= nchar(seq)] }
cut_positions <- function(seq, enz, circular = FALSE) { r <- enzymes[enzymes$name == enz, ]; st <- find_sites(seq, r$site, circular); st + r$cut - 1 }   # cut after this base
fragments <- function(seq, cuts, circular = FALSE) { L <- nchar(seq); cuts <- sort(unique(cuts))
  if (length(cuts) == 0) return(L)
  if (circular) { if (length(cuts) == 1) L else diff(c(cuts, cuts[1] + L)) } else diff(c(0, cuts, L)) }
digest <- function(seq, enz_vec, circular) { cuts <- unlist(lapply(enz_vec, function(e) cut_positions(seq, e, circular))); fragments(seq, cuts, circular) }

tm_basic <- function(p) { p <- toupper(p); n <- nchar(p); gc <- sum(strsplit(p, "")[[1]] %in% c("G", "C"))
  if (n < 14) 2 * (n - gc) + 4 * gc else 64.9 + 41 * (gc - 16.4) / n }          # Wallace rule below 14 nt, otherwise the basic GC formula
gc_pct <- function(p) 100 * mean(strsplit(toupper(p), "")[[1]] %in% c("G", "C"))
revcomp <- function(s) paste(rev(strsplit(chartr("ACGT", "TGCA", toupper(s)), "")[[1]]), collapse = "")
# Choose a primer length from the 3' side of the annealing region to reach a target Tm (default 60 C).
design_primer <- function(template_region, target_tm = 60, min_len = 18, max_len = 30) {
  for (n in min_len:max_len) { p <- substr(template_region, 1, n); if (tm_basic(p) >= target_tm) return(p) }; substr(template_region, 1, max_len) }

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 2) { insert <- read_fasta1(args[1]); plasmid <- read_fasta1(args[2]) } else {
    rnd <- function(n) paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = "")
    insert <- paste0("ATG", rnd(300), "GGATCC", rnd(300), "TAA")
    insert <- gsub("GAATTC|AAGCTT|CTCGAG|CATATG", "GCTGCT", insert)                         # keep the demo insert free of EcoRI, HindIII, XhoI, NdeI sites
    insert <- paste0(substr(insert, 1, 300), "GGATCC", substr(insert, 307, nchar(insert)))     # one internal BamHI site
    mcs <- "GAATTCGGATCCAAGCTTCTCGAGGTCGAC"; plasmid <- paste0(rnd(1500), mcs, rnd(1500)); plasmid <- gsub("GAATTC|GGATCC|AAGCTT|CTCGAG|GTCGAC", "GCTGCT", substr(plasmid, 1, 1500 - 1)) ; plasmid <- paste0(plasmid, "G", mcs, rnd(1500))
    cat("Practice data: synthetic insert and plasmid. The plasmid cloning site carries EcoRI, BamHI, HindIII, XhoI and SalI sites; some enzymes also cut elsewhere by chance\n\n") }
  cat(sprintf("Insert %d bp, plasmid %d bp (circular)\n\n", nchar(insert), nchar(plasmid)))
  tab <- data.frame(enzyme = enzymes$name, cuts_insert = sapply(enzymes$name, function(e) length(cut_positions(insert, e))), cuts_plasmid = sapply(enzymes$name, function(e) length(cut_positions(plasmid, e, TRUE))), row.names = NULL)
  print(tab, row.names = FALSE)
  good <- tab$enzyme[tab$cuts_plasmid == 1 & tab$cuts_insert == 0]; cat("\nEnzymes that cut the plasmid once and the insert not at all:", if (length(good)) good else "none", "\n")
  e1 <- good[1]; e2 <- good[2]; cat(sprintf("Chosen pair for directional cloning: %s and %s\n", e1, e2))
  cat(sprintf("\nPlasmid cut with %s and %s: fragment sizes %s bp\n", e1, e2, paste(sort(digest(plasmid, c(e1, e2), TRUE), decreasing = TRUE), collapse = ", ")))
  # Primers add a clamp and the chosen site to the 5' end of the annealing region.
  clamp <- "TTAGCA"; s1 <- enzymes$site[enzymes$name == e1]; s2 <- enzymes$site[enzymes$name == e2]
  fwd_anneal <- design_primer(insert); rev_anneal <- design_primer(revcomp(insert))
  fwd <- paste0(clamp, s1, fwd_anneal); rev <- paste0(clamp, s2, rev_anneal)
  cat(sprintf("\nForward primer 5'-%s-3'  (annealing part %d nt, Tm %.1f C, GC %.0f%%)\n", fwd, nchar(fwd_anneal), tm_basic(fwd_anneal), gc_pct(fwd_anneal)))
  cat(sprintf("Reverse primer 5'-%s-3'  (annealing part %d nt, Tm %.1f C, GC %.0f%%)\n", rev, nchar(rev_anneal), tm_basic(rev_anneal), gc_pct(rev_anneal)))
  cat(sprintf("Predicted PCR product: %d bp (insert %d bp plus %d bp of added ends)\n", nchar(insert) + 2 * (nchar(clamp)) + nchar(s1) + nchar(s2), nchar(insert), 2 * nchar(clamp) + nchar(s1) + nchar(s2)))
  cat("\nTm formula: Wallace rule (2 per A/T, 4 per G/C) under 14 nt; otherwise 64.9 + 41 x (G + C - 16.4) / length. Nearest-neighbour methods are more accurate.\n")
  cat("Insert digested with the enzymes that cut it (BamHI if present):", paste(digest(insert, "BamHI", FALSE), collapse = " + "), "bp\n")
  png("virtual_gel.png", 500, 600); sizes <- list(Plasmid_uncut = nchar(plasmid), Plasmid_digest = digest(plasmid, c(e1, e2), TRUE), Insert_BamHI = digest(insert, "BamHI", FALSE))
  plot(NA, xlim = c(0.5, 3.5), ylim = log10(c(100, 5000)), xaxt = "n", xlab = "", ylab = "log10 fragment size (bp)"); axis(1, 1:3, names(sizes)); for (i in 1:3) segments(i - 0.3, log10(sizes[[i]]), i + 0.3, log10(sizes[[i]]), lwd = 4); dev.off()
}
