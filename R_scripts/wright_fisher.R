# EVO-04  Drift and selection in simulated populations (Wright-Fisher model). Base R only: no packages to install.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
#
# One locus with two alleles, A and a, in a diploid population of N individuals (2N gene copies).
# Each generation: selection changes the frequency of A, then binomial sampling of 2N copies (drift).
#
# Run in RStudio: open this file and click Source. Run from a terminal: Rscript wright_fisher.R

simulate_wf <- function(N, p0 = 0.5, s = 0, h = 0.5, gens = 100, reps = 20, seed = 1) {
  # N     diploid population size (2N gene copies)
  # p0    starting frequency of allele A
  # s     selection coefficient favouring A (fitness AA = 1+s, Aa = 1+h*s, aa = 1)
  # h     dominance (0.5 = additive, 1 = A dominant, 0 = A recessive)
  # gens  number of generations
  # reps  number of replicate populations
  # Returns a matrix: rows = generations 0..gens, columns = replicates.
  set.seed(seed)
  p <- rep(p0, reps)
  traj <- matrix(NA_real_, nrow = gens + 1, ncol = reps)
  traj[1, ] <- p
  for (g in 1:gens) {
    q <- 1 - p
    w_AA <- 1 + s; w_Aa <- 1 + h * s; w_aa <- 1
    w_bar <- p^2 * w_AA + 2 * p * q * w_Aa + q^2 * w_aa
    p <- (p^2 * w_AA + p * q * w_Aa) / w_bar          # selection
    p <- rbinom(reps, 2 * N, p) / (2 * N)             # drift
    traj[g + 1, ] <- p
  }
  traj
}

summarise_wf <- function(traj) {
  final <- traj[nrow(traj), ]
  list(fixed = mean(final == 1),
       lost  = mean(final == 0),
       het   = rowMeans(2 * traj * (1 - traj)))       # mean 2pq for each generation
}

plot_wf <- function(traj, n_show = 20, main = "") {
  n_show <- min(n_show, ncol(traj))
  matplot(0:(nrow(traj) - 1), traj[, 1:n_show], type = "l", lty = 1,
          col = rgb(0.1, 0.35, 0.25, 0.6), ylim = c(0, 1),
          xlab = "Generation", ylab = "Frequency of allele A", main = main)
}

if (sys.nframe() == 0 || interactive()) {
  cat("Drift only, p0 = 0.5, 100 generations, 200 replicates\n")
  for (N in c(10, 100, 1000)) {
    traj <- simulate_wf(N, gens = 100, reps = 200, seed = 42)
    r <- summarise_wf(traj)
    theory <- 0.5 * (1 - 1 / (2 * N))^100
    cat(sprintf("N=%5d  fixed=%.2f  lost=%.2f  mean 2pq at gen 100 = %.3f  (theory %.3f)\n",
                N, r$fixed, r$lost, r$het[101], theory))
  }

  cat("\nNeutral allele, p0 = 0.05, N = 50, 2000 generations, 2000 replicates\n")
  traj <- simulate_wf(50, p0 = 0.05, gens = 2000, reps = 2000, seed = 3)
  r <- summarise_wf(traj)
  cat(sprintf("fixed=%.3f  lost=%.3f  (fixed should be close to 0.05)\n", r$fixed, r$lost))

  cat("\nSelection s = 0.05 (additive), p0 = 0.05, 600 generations, 200 replicates\n")
  for (N in c(20, 200, 2000)) {
    traj <- simulate_wf(N, p0 = 0.05, s = 0.05, gens = 600, reps = 200, seed = 7)
    r <- summarise_wf(traj)
    cat(sprintf("N=%5d  fixed=%.2f  lost=%.2f\n", N, r$fixed, r$lost))
  }
}
