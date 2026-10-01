# SYS-02  Barcode gap and species identification from aligned DNA sequences. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript sys02_barcode_gap.R [aligned.fasta]
# FASTA headers must read >Species_name|specimen_id (the part before the first | is the species).
# With no argument the script simulates 8 species with 6 specimens each along a tree, so it runs offline.
set.seed(658)

read_fasta <- function(file) { l <- readLines(file); h <- grep("^>", l); ends <- c(h[-1] - 1, length(l))
  setNames(mapply(function(s, e) toupper(paste(l[(s + 1):e], collapse = "")), h, ends), sub("^>", "", l[h])) }
# p-distance: proportion of differing sites, ignoring sites where either sequence has a gap or N.
p_dist_matrix <- function(seqs) { m <- do.call(rbind, strsplit(seqs, "")); n <- nrow(m); d <- matrix(0, n, n, dimnames = list(names(seqs), names(seqs)))
  for (i in 1:(n - 1)) for (j in (i + 1):n) { ok <- !(m[i, ] %in% c("-", "N") | m[j, ] %in% c("-", "N")); d[i, j] <- d[j, i] <- mean(m[i, ok] != m[j, ok]) }
  d }
simulate_barcodes <- function(n_species = 8, n_per = 6, len = 600, between = 0.10, within = 0.01) {
  root <- sample(c("A", "C", "G", "T"), len, TRUE); mutate <- function(s, rate) { k <- runif(len) < rate; s[k] <- sample(c("A", "C", "G", "T"), sum(k), TRUE); s }
  out <- character(0)
  for (sp in seq_len(n_species)) { sp_seq <- mutate(root, between); for (i in seq_len(n_per)) out[sprintf("Species%02d|spec%02d", sp, i)] <- paste(mutate(sp_seq, within), collapse = "") }
  out }
gap_summary <- function(d, sp) { same <- outer(sp, sp, "=="); diag(same) <- NA; w <- d[which(same, arr.ind = FALSE)]; b <- d[which(!same)]
  # per species: largest intraspecific distance vs smallest interspecific distance
  per <- t(sapply(unique(sp), function(s) { i <- sp == s; c(max_within = if (sum(i) > 1) max(d[i, i][upper.tri(d[i, i])]) else NA, min_between = min(d[i, !i])) }))
  list(within = w, between = b, per_species = per) }
identify_nn <- function(d, sp) sapply(seq_along(sp), function(i) { dd <- d[i, ]; dd[i] <- Inf; sp[which.min(dd)] })

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  seqs <- if (length(args) >= 1) read_fasta(args[1]) else { cat("Practice data: 8 simulated species x 6 specimens, 600 sites, about 1% within and 10% between species\n"); simulate_barcodes() }
  sp <- sub("\\|.*", "", names(seqs)); d <- p_dist_matrix(seqs)
  g <- gap_summary(d, sp)
  cat(sprintf("\nWithin-species distance: mean %.3f, maximum %.3f\nBetween-species distance: mean %.3f, minimum %.3f\n", mean(g$within), max(g$within), mean(g$between), min(g$between)))
  cat(sprintf("Barcode gap present (largest within < smallest between): %s; gap size %.3f\n", max(g$within) < min(g$between), min(g$between) - max(g$within)))
  cat("\nPer species: largest within-species distance and nearest other species\n"); print(round(g$per_species, 3))
  id <- identify_nn(d, sp); cat(sprintf("\nLeave-one-out nearest-neighbour identification: %.1f%% correct\n", 100 * mean(id == sp)))
  cat("\nWhen species are closer: identification accuracy and gap as between-species divergence shrinks (20 simulations each)\n")
  for (b in c(0.10, 0.05, 0.03, 0.02)) { r <- replicate(20, { s <- simulate_barcodes(between = b); spp <- sub("\\|.*", "", names(s)); dd <- p_dist_matrix(s); c(mean(identify_nn(dd, spp) == spp), max(gap_summary(dd, spp)$within) < min(gap_summary(dd, spp)$between)) })
    cat(sprintf("between-species divergence %.2f: accuracy %.2f, gap present in %.0f%% of runs\n", b, mean(r[1, ]), 100 * mean(r[2, ]))) }
  cat("\nAn unknown specimen (simulated from Species03) identified by nearest match: ")
  s3 <- names(seqs)[sp == "Species03"][1]; q <- seqs[[s3]]; refs <- seqs[names(seqs) != s3]; dq <- sapply(refs, function(r) mean(strsplit(q, "")[[1]] != strsplit(r, "")[[1]])); cat(sub("\\|.*", "", names(which.min(dq))), sprintf("(distance %.3f)\n", min(dq)))
}
