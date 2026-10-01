# GEN-06  Classifying mutations by their effect on the protein. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript gen06_classify_mutations.R [cds.fasta variants.csv]
# variants.csv columns: id, pos, ref, alt   (pos is 1-based within the coding sequence; for an insertion ref is "-",
# for a deletion alt is "-"; ref and alt may be several bases).
# With no arguments the script uses a SYNTHETIC coding sequence and variants so it runs offline.

# Standard genetic code, codons in TCAG order.
bases <- c("T", "C", "A", "G")
codons <- unlist(lapply(bases, function(a) unlist(lapply(bases, function(b) paste0(a, b, bases)))))
aa <- strsplit("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "")[[1]]
code <- setNames(aa, codons)

translate <- function(seq) { seq <- toupper(seq); n <- nchar(seq) %/% 3
  if (n == 0) return("")
  cod <- substring(seq, seq(1, by = 3, length.out = n), seq(3, by = 3, length.out = n))
  p <- unname(code[cod]); p[is.na(p)] <- "X"; stop_at <- which(p == "*")[1]
  paste(if (is.na(stop_at)) p else p[1:stop_at], collapse = "") }

apply_variant <- function(cds, pos, ref, alt) {
  ref <- if (ref == "-") "" else ref; alt <- if (alt == "-") "" else alt
  if (ref != "" && substr(cds, pos, pos + nchar(ref) - 1) != ref) stop("Reference base mismatch at position ", pos)
  if (ref == "") paste0(substr(cds, 1, pos), alt, substr(cds, pos + 1, nchar(cds)))   # insertion after pos
  else paste0(substr(cds, 1, pos - 1), alt, substr(cds, pos + nchar(ref), nchar(cds)))
}

classify <- function(cds, id, pos, ref, alt) {
  mut <- apply_variant(cds, pos, ref, alt); p0 <- translate(cds); p1 <- translate(mut)
  dl <- nchar(if (alt == "-") "" else alt) - nchar(if (ref == "-") "" else ref)
  s0 <- strsplit(p0, "")[[1]]; s1 <- strsplit(p1, "")[[1]]
  type <- if (dl %% 3 != 0) "frameshift"
          else if (dl != 0) "in-frame indel"
          else if (p0 == p1) "synonymous"
          else if (pos <= 3) "start loss"
          else if (!("*" %in% s1) && tail(s0, 1) == "*") "stop loss"
          else if ("*" %in% s1 && length(s1) < length(s0)) "nonsense"
          else "missense"
  cp <- (pos - 1) %/% 3 + 1
  chg <- if (type %in% c("missense", "nonsense", "synonymous", "start loss", "stop loss")) sprintf("p.%s%d%s", s0[cp], cp, s1[cp]) else ""
  data.frame(id = id, pos = pos, change = paste0(ref, ">", alt), effect = type, protein_change = chg,
             protein_length_before = length(s0), protein_length_after = length(s1))
}

# ---------- demo data (synthetic) ----------
demo_cds <- paste0("ATG", "GCTAAAGAAGTTCTGGGCGACATCCGTCACTGGAACCCAAGCTTAGTAGAAAGCGACCTG", "TACGGTAATGCAGTTATTGAGTAA")
demo_spec <- data.frame(id = c("synonymous_1", "missense_1", "missense_2", "nonsense_1", "start_loss", "deletion_1bp", "insertion_2bp", "deletion_3bp", "stop_loss"),
                        pos = c(6, 14, 20, 10, 2, 25, 30, 40, 87), len = c(1, 1, 1, 1, 1, 1, 0, 3, 1),
                        alt = c("C", "G", "T", "T", "C", "-", "GA", "-", "T"), stringsAsFactors = FALSE)
make_demo_vars <- function(cds, spec) { spec$ref <- ifelse(spec$len == 0, "-", substring(cds, spec$pos, spec$pos + pmax(spec$len, 1) - 1)); spec[, c("id", "pos", "ref", "alt")] }

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 2) { l <- readLines(args[1]); cds <- paste(l[!grepl("^>", l)], collapse = ""); vars <- read.csv(args[2], stringsAsFactors = FALSE)
  } else { cds <- demo_cds; vars <- make_demo_vars(cds, demo_spec); cat("Using the built-in SYNTHETIC coding sequence (", nchar(cds), " bases, ", nchar(translate(cds)) - 1, " amino acids)\n", sep = "") }
  cat("Reference protein:", translate(cds), "\n\n")
  options(width = 140)
  res <- do.call(rbind, lapply(seq_len(nrow(vars)), function(i) classify(cds, vars$id[i], vars$pos[i], vars$ref[i], vars$alt[i])))
  print(res, row.names = FALSE)
  cat("\nCounts by effect:\n"); print(table(res$effect))
  cat("\nCheck: genetic code translates ATGGCCTAA to", translate("ATGGCCTAA"), "and TGGTGATAG to", translate("TGGTGATAG"), "\n")
}
