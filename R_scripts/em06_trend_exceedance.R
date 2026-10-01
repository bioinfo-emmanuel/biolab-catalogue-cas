# EM-06  Trends and exceedances in public environmental monitoring data. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Input: CSV with columns date (YYYY-MM-DD) and value. With no file, uses R's built-in airquality
# data (daily ozone in ppb, New York, May to September 1973) so the script runs offline.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) { d <- read.csv(args[1]); d$date <- as.Date(d$date) } else {
  aq <- datasets::airquality; d <- data.frame(date = as.Date(sprintf("1973-%02d-%02d", aq$Month, aq$Day)), value = aq$Ozone) }
d <- d[!is.na(d$value), ]; d <- d[order(d$date), ]
guideline <- if (length(args) >= 2) as.numeric(args[2]) else 60   # set from a cited guideline for the same units

# Theil-Sen slope (median of pairwise slopes) and Mann-Kendall test (via Kendall's tau).
sen_slope <- function(t, y) { i <- combn(length(t), 2); s <- (y[i[2, ]] - y[i[1, ]]) / (t[i[2, ]] - t[i[1, ]]); median(s[is.finite(s)]) }
trend <- function(d) { t <- as.numeric(d$date); ct <- cor.test(t, d$value, method = "kendall", exact = FALSE)
  c(tau = unname(ct$estimate), p = ct$p.value, sen_per_day = sen_slope(t, d$value), ols_per_day = unname(coef(lm(value ~ t, d))[2])) }

# Exceedance counts by month, and the longest run of consecutive exceedance days.
exceedances <- function(d, guideline) { d$exc <- d$value > guideline; d$month <- format(d$date, "%Y-%m")
  list(by_month = aggregate(exc ~ month, d, function(x) c(days = length(x), exceed = sum(x))),
       longest_run = max(rle(d$exc)$lengths[rle(d$exc)$values], 0)) }

# Autocorrelation matters: lag-1 correlation of residuals from the trend line.
lag1 <- function(d) { r <- resid(lm(value ~ as.numeric(date), d)); cor(r[-1], r[-length(r)]) }

if (sys.nframe() == 0) {
  cat(sprintf("Records: %d  from %s to %s\nGuideline used: %g (same units as the data)\n", nrow(d), min(d$date), max(d$date), guideline))
  print(round(trend(d), 4))
  cat("\nMedian by month:\n"); print(round(tapply(d$value, format(d$date, "%m"), median), 1))
  e <- exceedances(d, guideline); cat("\nExceedance days by month\n"); print(e$by_month)
  cat("Longest run of consecutive exceedance days:", e$longest_run, "\n")
  cat(sprintf("\nLag-1 autocorrelation of residuals: %.2f (a value well above 0 means the trend p-value is too optimistic)\n", lag1(d)))
  cat("\nEffect of averaging: the same series as weekly means\n")
  d$wk <- format(d$date, "%Y-%U"); w <- aggregate(value ~ wk, d, mean); cat("weekly means above the guideline:", sum(w$value > guideline), "of", nrow(w), "weeks\n")
}
