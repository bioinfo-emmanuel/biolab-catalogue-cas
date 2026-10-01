# PAR-04  Mosquito larval survey and container indices. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript par04_larval_indices.R [containers.csv]
# containers.csv: area, house, container, container_type, has_water (1/0), positive (1 if larvae or pupae present), larvae_count (optional)
# Indices (World Health Organization definitions):
#   House index (HI)     = houses with at least one positive container / houses inspected x 100
#   Container index (CI) = positive containers / containers holding water x 100
#   Breteau index (BI)   = positive containers per 100 houses inspected
set.seed(4)
indices <- function(d) { houses <- length(unique(d$house)); w <- d[d$has_water == 1, ]; pos <- w[w$positive == 1, ]
  c(houses = houses, containers_with_water = nrow(w), positive_containers = nrow(pos), house_index = 100 * length(unique(pos$house)) / houses, container_index = 100 * nrow(pos) / nrow(w), breteau_index = 100 * nrow(pos) / houses) }
prop_ci <- function(x, n) 100 * binom.test(x, n)$conf.int
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) d <- read.csv(args[1]) else {
    types <- c("bucket", "drum", "tyre", "plant_saucer", "vase", "gutter", "coconut_shell", "can")
    mk <- function(area, nh, pos_rate) do.call(rbind, lapply(1:nh, function(h) { k <- rpois(1, 3) + 1; data.frame(area = area, house = paste0(area, h), container = 1:k, container_type = sample(types, k, TRUE), has_water = rbinom(k, 1, 0.6)) }))
    d <- rbind(mk("dorm_area", 40, 0.2), mk("canteen_area", 40, 0.35)); rate <- c(bucket = 0.25, drum = 0.30, tyre = 0.45, plant_saucer = 0.25, vase = 0.15, gutter = 0.20, coconut_shell = 0.40, can = 0.30); d$positive <- ifelse(d$has_water == 1 & runif(nrow(d)) < rate[d$container_type] * ifelse(d$area == "canteen_area", 1.3, 0.8), 1, 0); d$larvae_count <- ifelse(d$positive == 1, rpois(nrow(d), 12) + 1, 0)
    cat("Practice data: 40 invented premises in each of two invented areas; container types and positivity are simulated.\n\n") }
  res <- do.call(rbind, lapply(split(d, d$area), function(x) data.frame(area = x$area[1], t(round(indices(x), 1))))); print(res, row.names = FALSE)
  for (a in unique(d$area)) { x <- d[d$area == a, ]; w <- x[x$has_water == 1, ]; ci <- prop_ci(sum(w$positive), nrow(w)); pos_h <- length(unique(w$house[w$positive == 1])); ch <- prop_ci(pos_h, length(unique(x$house))); cat(sprintf("%s: container index 95%% CI %.0f to %.0f%%; house index 95%% CI %.0f to %.0f%%\n", a, ci[1], ci[2], ch[1], ch[2])) }
  ar <- unique(d$area); wtr <- d[d$has_water == 1, ]; f <- fisher.test(table(wtr$area, wtr$positive)); cat(sprintf("\nContainer positivity, %s vs %s: Fisher exact p = %.4f (odds ratio %.2f)\n", ar[1], ar[2], f$p.value, f$estimate))
  cat("\nContainer types ranked by the number of positive containers (all areas):\n"); tt <- aggregate(positive ~ container_type, wtr, function(v) c(with_water = length(v), positive = sum(v), pct = 100 * mean(v))); tt <- do.call(data.frame, tt); names(tt)[2:4] <- c("with_water", "positive", "pct_positive"); print(tt[order(-tt$positive), ], digits = 3, row.names = FALSE)
  cat("\nWhy several indices: house index ignores how many containers are positive; container index ignores the number of houses; Breteau index combines both but not larval density.\n")
  cat("Check: 8 positive houses of 40 gives HI =", 100 * 8 / 40, "; 10 positive of 50 wet containers gives CI =", 100 * 10 / 50, "; 10 positive containers in 40 houses gives BI =", 100 * 10 / 40, "\n")
}
