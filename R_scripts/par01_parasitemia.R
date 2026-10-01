# PAR-01  Estimating parasitemia from blood-film slides, with observer agreement. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript par01_parasitemia.R [counts.csv]
# counts.csv: slide, observer, method (thin or thick), parasitized (thin: parasitized red cells; thick: parasites), denominator (thin: red cells counted; thick: white cells counted)
# Thin film: parasitemia % = parasitized red cells / red cells counted x 100.
# Thick film: parasites per microlitre = parasites / white cells counted x 8000, using the conventional assumed count of 8000 white cells per microlitre.
set.seed(1)
thin_pct <- function(x, n) 100 * x / n
thick_per_ul <- function(p, w, wbc = 8000) p / w * wbc
wilson <- function(x, n) { b <- binom.test(x, n)$conf.int; 100 * b }
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) d <- read.csv(args[1]) else {
    truth <- c(0.4, 1.2, 3.0, 6.0, 12.0) / 100; rows <- list()
    for (i in seq_along(truth)) for (ob in c("A", "B")) { n <- 500; rows[[length(rows) + 1]] <- data.frame(slide = paste0("S", i), observer = ob, method = "thin", parasitized = rbinom(1, n, truth[i] * ifelse(ob == "B", 1.1, 1)), denominator = n) }
    dens <- c(100, 500, 2000, 8000, 20000)       # invented true parasite densities per microlitre for the thick films
    for (i in seq_along(dens)) for (ob in c("A", "B")) { w <- 200; rows[[length(rows) + 1]] <- data.frame(slide = paste0("T", i), observer = ob, method = "thick", parasitized = rpois(1, dens[i] / 8000 * w), denominator = w) }
    d <- do.call(rbind, rows); cat("Practice data: 5 invented thin films (true parasitemia 0.4 to 12 percent) and 5 invented thick films (100 to 20000 parasites per microlitre), each counted by two observers: 500 red cells per thin film and 200 white cells per thick film. The counts are simulated, not from real slides.\n\n") }
  th <- d[d$method == "thin", ]; th$pct <- thin_pct(th$parasitized, th$denominator); ci <- t(mapply(wilson, th$parasitized, th$denominator)); th$lo <- ci[, 1]; th$hi <- ci[, 2]; print(format(th[, c("slide", "observer", "parasitized", "denominator", "pct", "lo", "hi")], digits = 3), row.names = FALSE)
  w <- reshape(th[, c("slide", "observer", "pct")], idvar = "slide", timevar = "observer", direction = "wide"); df <- w[[3]] - w[[2]]
  cat(sprintf("\nObserver agreement on thin-film parasitemia: mean difference %.2f percentage points (B minus A), limits of agreement %.2f to %.2f; correlation %.3f\n", mean(df), mean(df) - 1.96 * sd(df), mean(df) + 1.96 * sd(df), cor(w[[2]], w[[3]])))
  cat("\nHow many red cells to count? Width of the 95% interval at 2% true parasitemia:\n"); for (n in c(100, 200, 500, 1000, 2000)) { x <- round(0.02 * n); ci <- wilson(x, n); cat(sprintf("n = %4d cells: %.2f%% (%.2f to %.2f), relative width %.0f%%\n", n, 100 * x / n, ci[1], ci[2], 100 * (ci[2] - ci[1]) / (100 * x / n))) }
  tk <- d[d$method == "thick", ]; if (nrow(tk)) { tk$per_ul <- thick_per_ul(tk$parasitized, tk$denominator); cat("\nThick film (parasites per microlitre, using 8000 white cells per microlitre), mean of the two observers:\n"); ag <- aggregate(per_ul ~ slide, tk, mean); ag$per_ul <- round(ag$per_ul, 0); print(ag, row.names = FALSE) }
  cat("\nCheck: 12 parasitized cells among 500 red cells is", round(thin_pct(12, 500), 1), "percent; 30 parasites against 200 white cells is", thick_per_ul(30, 200), "per microlitre\n")
}
