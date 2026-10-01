# HG-01  Genetic counseling: Bayesian carrier probability. Base R only.
# Compiled by Emmanuel G. Lozano, Department of Biology, College of Arts and Sciences, Pampanga State University. CC BY-NC 4.0: credit required, non-commercial use only.
# Prior, conditional, joint, posterior for an unaffected person, using offspring outcomes.

posterior <- function(prior_carrier, p_data_if_carrier, p_data_if_not) {
  j1 <- prior_carrier * p_data_if_carrier; j0 <- (1 - prior_carrier) * p_data_if_not
  c(prior = prior_carrier, joint_carrier = j1, joint_not = j0, posterior = j1 / (j1 + j0))
}

# Autosomal recessive. Unaffected sibling of an affected person (parents both carriers): prior 2/3.
# With a partner of carrier probability q, the child is affected with probability carrier * q * 1/4.
ar_sibling_risk <- function(q_partner, n_unaffected_children = 0, prior = 2 / 3) {
  # If she already has n unaffected children with a partner who is a carrier with probability q:
  # P(child unaffected | she is carrier) = 1 - q/4 per child; if she is not a carrier, 1.
  post <- posterior(prior, (1 - q_partner / 4)^n_unaffected_children, 1)
  carrier <- post["posterior"]
  c(carrier = unname(carrier), risk_next_child = unname(carrier * q_partner / 4))
}

# X-linked recessive. A woman whose mother is a known carrier has prior 1/2.
# Each unaffected son has probability 1/2 if she is a carrier and 1 if she is not.
xlr_daughter <- function(n_unaffected_sons, prior = 1 / 2) posterior(prior, 0.5^n_unaffected_sons, 1)

# Carrier frequency from disease incidence q^2 under Hardy-Weinberg: carriers = 2pq.
carrier_freq <- function(incidence) { q <- sqrt(incidence); 2 * q * (1 - q) }

if (sys.nframe() == 0) {
  cat("Autosomal recessive: healthy sister of an affected person, partner carrier probability 1/25\n")
  for (n in 0:3) { r <- ar_sibling_risk(1 / 25, n); cat(sprintf("%d unaffected children: P(she is a carrier) = %.3f, risk for next child = %.4f\n", n, r["carrier"], r["risk_next_child"])) }
  cat("\nSame sister, but her partner is a known carrier (q = 1): unaffected children now lower her carrier probability\n")
  for (n in 0:3) { r <- ar_sibling_risk(1, n); cat(sprintf("%d unaffected children: P(she is a carrier) = %.3f, risk for next child = %.4f\n", n, r["carrier"], r["risk_next_child"])) }
  cat("\nX-linked recessive: mother is a known carrier; the consultand is her daughter\n")
  for (k in 0:4) { r <- xlr_daughter(k); cat(sprintf("%d unaffected sons: P(carrier) = %.3f\n", k, r["posterior"])) }
  cat("\nCarrier frequency from incidence (Hardy-Weinberg)\n")
  for (inc in c(1e-4, 4e-4, 1e-3)) cat(sprintf("incidence 1 in %.0f -> carrier frequency about 1 in %.0f\n", 1 / inc, 1 / carrier_freq(inc)))
  cat("\nTable of the joint calculation for 2 unaffected sons, X-linked:\n"); print(round(xlr_daughter(2), 3))
}
