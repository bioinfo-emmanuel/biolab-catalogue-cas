# EVO-02  Morphometrics and selection in a wild population. Base R only: no packages to install.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
#
# YOUR DATA: a CSV with one row per measurement and these columns
#   specimen, site, measure, height, width
#   optional: aperture (any further trait), status (live or dead; dead = predated shells)
#   specimen    a unique id for each specimen
#   site        a name for each site (two sites are compared)
#   measure     1 for the first measurement of the specimen, 2 for the repeat by another person on another day
#   height, width  in mm, or any unit used consistently
# The specimen is the replicate. Measure every specimen twice.
#
# Run in RStudio: set data_file below and click Source. From a terminal: Rscript morphometrics_analysis.R mydata.csv
# With no file the script makes PRACTICE data (simulated). Never report them as results.

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1) args[1] else NA

if (is.na(data_file)) {
  message("No data file given: using simulated PRACTICE data.")
  set.seed(12)
  mk <- function(site, n, h_mean, status = "live") {
    h <- rnorm(n, h_mean, 2.2); w <- 0.78 * h^1.0 * exp(rnorm(n, 0, 0.04)); a <- 0.60 * h * exp(rnorm(n, 0, 0.05))
    id <- paste0(site, "_", status, "_", seq_len(n))
    rbind(data.frame(specimen = id, site = site, measure = 1, height = round(h + rnorm(n, 0, 0.3), 1), width = round(w + rnorm(n, 0, 0.3), 1), aperture = round(a + rnorm(n, 0, 0.3), 1), status = status),
          data.frame(specimen = id, site = site, measure = 2, height = round(h + rnorm(n, 0, 0.3), 1), width = round(w + rnorm(n, 0, 0.3), 1), aperture = round(a + rnorm(n, 0, 0.3), 1), status = status))
  }
  d <- rbind(mk("Site_A", 30, 22), mk("Site_B", 30, 25), mk("Site_B", 15, 27.5, "dead"))
} else {
  d <- read.csv(data_file, stringsAsFactors = FALSE)
}
traits <- intersect(c("height", "width", "aperture"), names(d))
if (!"status" %in% names(d)) d$status <- "live"

# ---- 1. repeatability of each trait (two measurements per specimen) --------------------
cat("\n1. Repeatability r = (MSamong - MSwithin) / (MSamong + MSwithin), from a one-way ANOVA with specimens as groups\n")
rep_tab <- sapply(traits, function(tr) { a <- anova(aov(d[[tr]] ~ factor(d$specimen))); msa <- a[1, "Mean Sq"]; msw <- a[2, "Mean Sq"]; (msa - msw) / (msa + msw) })
print(round(rep_tab, 3))
cat("Values below about 0.9 mean the measurement error is large relative to real differences.\n")

# ---- 2. one value per specimen: the mean of its two measurements --------------------------
s <- aggregate(d[, traits, drop = FALSE], by = list(specimen = d$specimen, site = d$site, status = d$status), FUN = mean)
live <- s[s$status == "live", ]

cat("\n2. Mean, sd, and coefficient of variation by site (live specimens)\n")
for (tr in traits) print(do.call(rbind, lapply(split(live, live$site), function(x) data.frame(trait = tr, site = x$site[1], n = nrow(x),
    mean = round(mean(x[[tr]]), 2), sd = round(sd(x[[tr]]), 2), cv_pct = round(100 * sd(x[[tr]]) / mean(x[[tr]]), 1)))), row.names = FALSE)

cat("\n3. Two-site comparison (Welch t-test, live specimens)\n")
sites <- unique(live$site)
if (length(sites) == 2) for (tr in traits) { t <- t.test(live[[tr]][live$site == sites[1]], live[[tr]][live$site == sites[2]])
  cat(sprintf("%-9s difference (%s minus %s) = %.2f, 95%% CI %.2f to %.2f, p = %.2e\n", tr, sites[1], sites[2], unname(t$estimate[1] - t$estimate[2]), t$conf.int[1], t$conf.int[2], t$p.value)) }

cat("\n4. Allometry: slope of log10 width on log10 height (a slope whose 95% interval includes 1 is consistent with isometry)\n")
if (all(c("height", "width") %in% traits)) { f <- lm(log10(width) ~ log10(height), data = live); ci <- confint(f)[2, ]
  cat(sprintf("slope = %.2f (95%% CI %.2f to %.2f), R squared %.2f\n", coef(f)[2], ci[1], ci[2], summary(f)$r.squared)) }

cat("\n5. Selection differential (optional; needs live and dead shells from the same site)\n")
if (any(s$status == "dead")) { for (site in unique(s$site[s$status == "dead"])) { x <- s[s$site == site, ]; tr <- traits[1]
    surv <- x[[tr]][x$status == "live"]; coh <- x[[tr]]; S <- mean(surv) - mean(coh)
    cat(sprintf("%s, %s: S = %.2f, standardized i = S / sd(cohort) = %.2f (live n = %d, cohort n = %d)\n", site, tr, S, S / sd(coh), length(surv), length(coh))) }
  cat("This assumes dead shells are the individuals removed by selection and nothing else, and that trait size did not change after death.\n") } else cat("No dead specimens in the data.\n")

par(mfrow = c(1, 2)); boxplot(height ~ site, data = live, ylab = "Height", main = "Height by site")
plot(live$height, live$width, log = "xy", pch = 19, xlab = "Height", ylab = "Width", main = "Allometry (log axes)")
