# ECO-02  Mark-recapture estimates and a test of the assumptions. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Usage: Rscript eco02_mark_recapture.R [M C R]      (M marked in the first sample, C caught in the second, R recaptured marked)
# With no arguments the script runs a simulation of a bucket of beans (or a real population) with a known size N.
set.seed(2)
lincoln_petersen <- function(M, C, R) M * C / R
chapman <- function(M, C, R) (M + 1) * (C + 1) / (R + 1) - 1
chapman_var <- function(M, C, R) (M + 1) * (C + 1) * (M - R) * (C - R) / ((R + 1)^2 * (R + 2))
chapman_ci <- function(M, C, R, level = 0.95) { est <- chapman(M, C, R); se <- sqrt(chapman_var(M, C, R)); z <- qnorm(1 - (1 - level) / 2); c(estimate = est, lower = est - z * se, upper = est + z * se) }

# One simulated two-visit study. Options break the assumptions:
#  trap_response > 1 : marked animals are more likely (trap-happy) or < 1 less likely (trap-shy) to be caught again
#  mark_loss         : fraction of marks lost between visits
#  deaths, births    : fraction of the population that dies or is added between visits (population not closed)
simulate_study <- function(N = 500, M = 100, C = 100, trap_response = 1, mark_loss = 0, deaths = 0, births = 0) {
  marked <- rep(FALSE, N); marked[sample.int(N, M)] <- TRUE
  alive <- rep(TRUE, N); dead <- sample.int(N, round(deaths * N)); alive[dead] <- FALSE
  pop <- c(marked & alive, rep(FALSE, round(births * N))); alive2 <- c(alive, rep(TRUE, round(births * N)))
  keep <- which(alive2); w <- ifelse(pop[keep], trap_response, 1); take <- keep[sample.int(length(keep), min(C, length(keep)), prob = w)]
  R <- sum(pop[take] & runif(length(take)) > mark_loss); c(M = M, C = length(take), R = R) }

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 3) { M <- as.numeric(args[1]); C <- as.numeric(args[2]); R <- as.numeric(args[3])
    cat(sprintf("Lincoln-Petersen N = %.1f\nChapman N = %.1f, 95%% interval %.1f to %.1f\n", lincoln_petersen(M, C, R), chapman(M, C, R), chapman_ci(M, C, R)[2], chapman_ci(M, C, R)[3])) } else {
    cat("Worked example: 60 marked, 50 caught later, 12 of them marked\n"); print(round(c(lincoln_petersen = lincoln_petersen(60, 50, 12), chapman_ci(60, 50, 12)), 1))
    cat("\nBucket of beans (true N = 500): 100 marked, 100 drawn later. 1000 simulated studies for each condition; mean of the Chapman estimate and share of 95% intervals that contain 500\n")
    run <- function(label, ...) { opts <- list(...); est <- t(replicate(1000, { s <- do.call(simulate_study, opts); if (s["R"] < 2) c(NA, NA, NA) else chapman_ci(s["M"], s["C"], s["R"]) })); est <- est[!is.na(est[, 1]), ]
      cat(sprintf("%-38s mean estimate %6.1f, coverage %.2f\n", label, mean(est[, 1]), mean(est[, 2] <= 500 & est[, 3] >= 500))) }
    run("all assumptions met"); run("trap-happy (marked 3 x as catchable)", trap_response = 3); run("trap-shy (marked 0.3 x as catchable)", trap_response = 0.3)
    run("20% of marks lost", mark_loss = 0.2); run("10% of animals die between visits", deaths = 0.10); run("20% of population added by births", births = 0.20)
    cat("\nEffect of sample size (assumptions met): 1000 studies, mean estimate and share of intervals containing N\n")
    for (mc in c(20, 50, 100, 200)) run(sprintf("M = C = %d", mc), M = mc, C = mc)
    cat("\nWith M = C = 20 and N = 500 the expected number of recaptures is only", 20 * 20 / 500, "so many studies have R = 0 or 1 and cannot be analysed.\n") }
}
