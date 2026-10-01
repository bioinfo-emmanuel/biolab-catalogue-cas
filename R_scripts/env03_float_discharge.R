# ENV-03  Stream velocity and discharge by the float method. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript env03_float_discharge.R [times.csv depths.csv reach_m width_m [k]]
# times.csv: float, time_s   (travel time of each float over the reach)   depths.csv: station_m, depth_m   (depth at stations across the channel, including 0 at each bank)
# k is the coefficient that converts surface velocity to mean velocity. Default 0.85 (a commonly used value; guidance ranges from about 0.8 for rough beds to 0.9 for smooth beds, and 0.66 to 0.75 has been quoted for shallow channels).
set.seed(77)
area_trapezoid <- function(x, z) sum(diff(x) * (z[-1] + z[-length(z)]) / 2)
discharge <- function(times, reach_m, area, k = 0.85) { vs <- reach_m / times; c(v_surface = mean(vs), v_mean = k * mean(vs), Q = k * mean(vs) * area) }

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE); k <- if (length(args) >= 5) as.numeric(args[5]) else 0.85
  if (length(args) >= 4) { tm <- read.csv(args[1]); dp <- read.csv(args[2]); reach <- as.numeric(args[3]) } else {
    reach <- 20; tm <- data.frame(float = 1:6, time_s = round(c(28.4, 25.1, 27.9, 24.6, 26.3, 29.0) + rnorm(6, 0, 0.3), 1))
    dp <- data.frame(station_m = seq(0, 6, by = 1), depth_m = c(0, 0.22, 0.38, 0.45, 0.41, 0.26, 0)); cat("Practice data: a 20 m reach, 6 floats, a 6 m wide channel measured every 1 m (invented)\n\n") }
  A <- area_trapezoid(dp$station_m, dp$depth_m); r <- discharge(tm$time_s, reach, A, k)
  cat(sprintf("Cross-section area %.2f m2 (trapezoid rule over %d stations); mean depth %.2f m\n", A, nrow(dp), A / diff(range(dp$station_m))))
  cat(sprintf("Surface velocity %.3f m/s (sd of individual floats %.3f); mean velocity with k = %.2f: %.3f m/s\nDischarge Q = %.3f m3/s = %.0f L/s\n", r["v_surface"], sd(reach / tm$time_s), k, r["v_mean"], r["Q"], 1000 * r["Q"]))
  cat(sprintf("Shortest float time %.1f s; guidance is a reach that takes floats more than about 20 s when possible\n", min(tm$time_s)))
  cat("\nSensitivity of discharge to k:\n"); for (kk in c(0.66, 0.75, 0.85, 0.9)) cat(sprintf("k = %.2f: Q = %.3f m3/s\n", kk, discharge(tm$time_s, reach, A, kk)["Q"]))
  cat("\nUncertainty from timing and area (bootstrap over floats and 5 percent random error in each depth):\n")
  q <- replicate(2000, { t <- sample(tm$time_s, replace = TRUE); z <- dp$depth_m * rnorm(nrow(dp), 1, 0.05); discharge(t, reach, area_trapezoid(dp$station_m, z), k)["Q"] }); cat(sprintf("Q median %.3f m3/s, 95%% interval %.3f to %.3f\n", median(q), quantile(q, .025), quantile(q, .975)))
}
