# BCH-06  Backbone dihedral angles and a Ramachandran plot from a PDB file. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript bch06_ramachandran.R [structure.pdb]
# Reads ATOM records of the first model (first alternate location only) and computes phi and psi for each residue.
# With no argument the script BUILDS two ideal test chains (an alpha helix and a beta strand) with known angles,
# checks that the calculation recovers them, and uses them as the practice structure, so it runs offline.

read_pdb_backbone <- function(lines) {
  a <- lines[grepl("^ATOM", lines)]
  end <- which(grepl("^ENDMDL", lines))[1]; if (!is.na(end)) a <- lines[seq_len(end - 1)][grepl("^ATOM", lines[seq_len(end - 1)])]
  d <- data.frame(name = trimws(substr(a, 13, 16)), alt = substr(a, 17, 17), res = trimws(substr(a, 18, 20)), chain = substr(a, 22, 22),
                  num = as.integer(substr(a, 23, 26)), x = as.numeric(substr(a, 31, 38)), y = as.numeric(substr(a, 39, 46)), z = as.numeric(substr(a, 47, 54)), stringsAsFactors = FALSE)
  d <- d[d$name %in% c("N", "CA", "C") & d$alt %in% c(" ", "A"), ]; d }

cross <- function(a, b) c(a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1])
dot <- function(a, b) sum(a * b)
dihedral <- function(p1, p2, p3, p4) { b1 <- p2 - p1; b2 <- p3 - p2; b3 <- p4 - p3; n1 <- cross(b1, b2); n2 <- cross(b2, b3)
  atan2(dot(cross(n1, n2), b2 / sqrt(dot(b2, b2))), dot(n1, n2)) * 180 / pi }

phi_psi <- function(bb) {
  out <- list()
  for (ch in unique(bb$chain)) { s <- bb[bb$chain == ch, ]; nums <- sort(unique(s$num))
    get <- function(n, at) { r <- s[s$num == n & s$name == at, ]; if (nrow(r) == 1) c(r$x, r$y, r$z) else NULL }
    for (n in nums) { N <- get(n, "N"); CA <- get(n, "CA"); C <- get(n, "C"); Cp <- get(n - 1, "C"); Nn <- get(n + 1, "N")
      phi <- if (!is.null(Cp) && !is.null(N) && !is.null(CA) && !is.null(C)) dihedral(Cp, N, CA, C) else NA
      psi <- if (!is.null(Nn) && !is.null(N) && !is.null(CA) && !is.null(C)) dihedral(N, CA, C, Nn) else NA
      out[[length(out) + 1]] <- data.frame(chain = ch, num = n, res = s$res[s$num == n][1], phi = phi, psi = psi) } }
  do.call(rbind, out) }

# Coarse teaching regions (rectangles, degrees). These are simplified boxes, not the validated contours used by structure-checking software.
region <- function(phi, psi) ifelse(is.na(phi) | is.na(psi), NA,
  ifelse(phi < 0 & phi > -160 & psi > -100 & psi < 50, "helix",
  ifelse(phi < -45 & psi > 90 | (phi < -45 & psi < -170), "strand",
  ifelse(phi > 0 & phi < 120 & psi > -30 & psi < 100, "left-handed", "other"))))

# ---- build ideal chains by placing each atom from the previous three (NeRF) ----
place <- function(A, B, C, r, theta_deg, phi_deg) {
  th <- theta_deg * pi / 180; ph <- phi_deg * pi / 180
  bc <- (C - B) / sqrt(dot(C - B, C - B)); n <- cross(B - A, bc); n <- n / sqrt(dot(n, n)); m <- cross(n, bc)
  d2 <- c(-r * cos(th), r * sin(th) * cos(ph), r * sin(th) * sin(ph)); C + d2[1] * bc + d2[2] * m + d2[3] * n }
build_chain <- function(n, phi, psi, chain = "A") {
  # ideal geometry (Engh and Huber style values): N-CA 1.458, CA-C 1.525, C-N 1.329 A; angles N-CA-C 111.2, CA-C-N 116.2, C-N-CA 121.7
  P <- list(N = c(0, 0, 0), CA = c(1.458, 0, 0)); ang <- 111.2 * pi / 180
  P$C <- P$CA + 1.525 * c(-cos(ang), sin(ang), 0)
  atoms <- list(N = P$N, CA = P$CA, C = P$C); coords <- list(); add <- function(i, nm, p) coords[[length(coords) + 1]] <<- list(i = i, nm = nm, p = p)
  add(1, "N", P$N); add(1, "CA", P$CA); add(1, "C", P$C); pN <- P$N; pCA <- P$CA; pC <- P$C
  for (i in 2:n) {
    Nn <- place(pN, pCA, pC, 1.329, 116.2, psi)               # psi of residue i-1 sets N(i)
    CAn <- place(pCA, pC, Nn, 1.458, 121.7, 180)              # omega = 180 (trans)
    Cn <- place(pC, Nn, CAn, 1.525, 111.2, phi)               # phi of residue i sets C(i)
    add(i, "N", Nn); add(i, "CA", CAn); add(i, "C", Cn); pN <- Nn; pCA <- CAn; pC <- Cn }
  sprintf("ATOM  %5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f  1.00  0.00           %s", seq_along(coords), sapply(coords, function(a) paste0(" ", a$nm)), "ALA", chain,
          sapply(coords, `[[`, "i"), sapply(coords, function(a) a$p[1]), sapply(coords, function(a) a$p[2]), sapply(coords, function(a) a$p[3]), substr(sapply(coords, `[[`, "nm"), 1, 1)) }

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  cat("Check of the dihedral function on a hand-built case: p1 = (1,0,0), p2 = (0,0,0), p3 = (0,0,1), p4 = (cos60, sin60, 1)\n")
  cat("expected +60, got", round(dihedral(c(1, 0, 0), c(0, 0, 0), c(0, 0, 1), c(cos(pi / 3), sin(pi / 3), 1)), 3), "\n")
  cat("expected 180 (trans), got", round(dihedral(c(1, 0, 0), c(0, 0, 0), c(0, 0, 1), c(-1, 0, 1)), 3), "\n\n")
  if (length(args) >= 1) { pdb <- readLines(args[1]) } else {
    pdb <- c(build_chain(20, -57, -47, "A"), build_chain(20, -120, 130, "B")); cat("Practice structure: two built chains, a 20-residue helix (phi -57, psi -47) and a 20-residue strand (phi -120, psi 130)\n") }
  bb <- read_pdb_backbone(pdb); a <- phi_psi(bb); a$region <- region(a$phi, a$psi)
  if (length(args) < 1) { h <- a[a$chain == "A" & !is.na(a$phi) & !is.na(a$psi), ]; s <- a[a$chain == "B" & !is.na(a$phi) & !is.na(a$psi), ]
    cat(sprintf("Recovered helix angles: phi %.1f (sd %.2f), psi %.1f (sd %.2f)\nRecovered strand angles: phi %.1f (sd %.2f), psi %.1f (sd %.2f)\n", mean(h$phi), sd(h$phi), mean(h$psi), sd(h$psi), mean(s$phi), sd(s$phi), mean(s$psi), sd(s$psi))) }
  cat(sprintf("\nResidues with both angles: %d of %d\n", sum(!is.na(a$phi) & !is.na(a$psi)), nrow(a)))
  cat("Residues per region:\n"); print(table(a$region, useNA = "no"))
  cat("\nFirst rows:\n"); print(head(a, 5), digits = 4, row.names = FALSE)
  png("ramachandran.png", 600, 600); plot(a$phi, a$psi, xlim = c(-180, 180), ylim = c(-180, 180), xlab = "phi", ylab = "psi", pch = 19, col = as.integer(factor(a$chain)) + 1); abline(h = 0, v = 0, lty = 3); dev.off()
}
