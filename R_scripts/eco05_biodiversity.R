# ECO-05  Comparing biodiversity from occurrence records. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript eco05_biodiversity.R [records.csv]
# records.csv columns: site (or habitat/region), species, and optionally count (default 1 record per row).
# With no argument the script simulates two communities with the same richness but different evenness, and
# a sampling-effort bias, so it runs offline.
set.seed(215)

abund_table <- function(d) { if (!"count" %in% names(d)) d$count <- 1; tapply(d$count, list(d$site, d$species), sum, default = 0) }
shannon <- function(x) { p <- x[x > 0] / sum(x); -sum(p * log(p)) }
simpson <- function(x) { p <- x / sum(x); 1 - sum(p^2) }                     # Gini-Simpson
chao1 <- function(x) { x <- x[x > 0]; f1 <- sum(x == 1); f2 <- sum(x == 2); S <- length(x)
  if (f2 > 0) S + f1^2 / (2 * f2) else S + f1 * (f1 - 1) / 2 }              # bias-corrected when f2 = 0

# Individual-based rarefaction (Hurlbert 1971): expected species in a random sample of m individuals.
rarefy <- function(x, m) { x <- x[x > 0]; N <- sum(x)
  if (m > N) return(NA)
  sum(1 - exp(lchoose(N - x, m) - lchoose(N, m))) }
accumulation <- function(x, reps = 100) { ind <- rep(seq_along(x), x); N <- length(ind)
  rowMeans(replicate(reps, { s <- sample(ind); cumsum(!duplicated(s)) })) }

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) dat <- read.csv(args[1]) else {
    sp <- paste0("sp", 1:40)
    even <- rmultinom(1, 400, rep(1, 40))[, 1]                               # even community
    uneven <- rmultinom(1, 400, 0.85^(0:39) / sum(0.85^(0:39)))[, 1]        # dominated community
    small_sample <- rmultinom(1, 60, rep(1, 40))[, 1]                        # same community as 'even', fewer records
    mk <- function(site, cnt) data.frame(site = site, species = sp[cnt > 0], count = cnt[cnt > 0])
    dat <- rbind(mk("Habitat_even", even), mk("Habitat_uneven", uneven), mk("Habitat_even_lowEffort", small_sample))
    cat("Practice data: 40 possible species; two habitats with 400 records, one habitat sampled with only 60 records\n\n") }
  tab <- abund_table(dat)
  res <- data.frame(site = rownames(tab), records = rowSums(tab), observed_S = rowSums(tab > 0),
                    shannon = apply(tab, 1, shannon), simpson = apply(tab, 1, simpson), chao1 = apply(tab, 1, chao1), singletons = rowSums(tab == 1))
  print(format(res, digits = 3), row.names = FALSE)
  m <- min(rowSums(tab)); cat(sprintf("\nRarefied richness at the smallest sample size (%d records):\n", m))
  print(round(setNames(apply(tab, 1, rarefy, m = m), rownames(tab)), 1))
  cat("\nExpected species at 20, 40 and 60 records:\n"); print(round(t(sapply(c(20, 40, 60), function(k) apply(tab, 1, rarefy, m = k))), 1))
  cat("\nBootstrap 95% interval for Shannon diversity of each site (500 resamples of records)\n")
  for (s in rownames(tab)) { ind <- rep(colnames(tab), tab[s, ]); b <- replicate(500, shannon(table(sample(ind, replace = TRUE))))
    cat(sprintf("%-24s %.2f to %.2f\n", s, quantile(b, .025), quantile(b, .975))) }
}
