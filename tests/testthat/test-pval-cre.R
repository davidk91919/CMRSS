# Test comparing CMRSS pval_comb_block against RIQITE::pval_quantile.
#
# Requires the RIQITE package (install from GitHub: li-xinran/RIQITE).
# Skipped if RIQITE or a solver is unavailable.
#
# CMRSS tests H_{k,c}^treat (treated-only) with k in 1..sum(Z); RIQITE tests
# H_{k,c}     (all-units)    with k in 1..n.  These are different hypotheses,
# but their LP minima for the rank-sum statistic coincide under the translation
#
#     k_R (RIQITE, all-units) = k_C (CMRSS, treated-only) + (n - m).
#
# Reason: in the LP, xi[i] only enters the statistic through Z[i] * xi[i], so
# controls' xi is unconstrained.  Assigning all (n - m) controls to xi = c
# absorbs (n - m) of RIQITE's all-units constraint slack, leaving exactly
# k_C = k_R - (n - m) treated units constrained.  With a shared Z.perm matrix
# the two packages return identical p-values, except when a null draw
# ties the observed value (see the Wilcoxon case below).

test_that("pval_comb_block matches RIQITE::pval_quantile under shared null and k translation", {
  # Skip if RIQITE is not available
  skip_if_not(
    requireNamespace("RIQITE", quietly = TRUE),
    "RIQITE package not available (install from GitHub: li-xinran/RIQITE)"
  )

  # Also need a solver for CMRSS
  skip_if_not(
    solver_available("highs") || solver_available("gurobi"),
    "Neither HiGHS nor Gurobi solver is available"
  )

  set.seed(54321) # Set seed for reproducibility

  n <- 50  # Fixed size for reproducibility
  m <- 25  # Number of treated units
  Z <- c(rep(1, m), rep(0, n - m))
  Z <- sample(Z)  # Shuffle

  Y <- rnorm(n)
  k_C   <- floor(0.8 * m)        # treated-only k for CMRSS, = 20
  k_R   <- k_C + (n - m)         # equivalent all-units k for RIQITE, = 45
  c_val <- 0

  r <- 3  # Stephenson parameter
  block <- factor(rep(1, n))

  # Shared permutation matrix (n x nperm) used by both packages so the empirical
  # null distributions match cell-for-cell.
  nperm <- 1000
  set.seed(98765)
  Z.perm <- RIQITE::assign_CRE(n, m, nperm)

  method.list.riqite.wil <- list(name = "Wilcoxon")
  method.list.riqite.ste <- list(name = "Stephenson", s = r)

  methods.list.cmrss.wil <- list(list(list(name = "Wilcoxon", scale = FALSE)))
  methods.list.cmrss.ste <- list(list(list(name = "Stephenson", scale = FALSE, s = r)))

  block.sum <- summary_block(Z, block)
  weight    <- weight_scheme(block.sum, "asymp.opt")

  # --- Wilcoxon ---
  riq_null_wil <- RIQITE::null_dist(n, m, nperm = nperm, Z.perm = Z.perm,
                                    method.list = method.list.riqite.wil)
  scores_wil   <- list(score_all_blocks(block.sum$nb, methods.list.cmrss.wil[[1]]))
  cms_null_wil <- com_null_dist_block(Z, block, methods.list.cmrss.wil, scores_wil,
                                      null.max = nperm, weight, block.sum,
                                      Z.perm = Z.perm)

  riqite_wil <- RIQITE::pval_quantile(
    Z = Z, Y = Y, k = k_R, c = c_val,
    method.list = method.list.riqite.wil,
    stat.null = riq_null_wil, switch = FALSE
  )

  cmrss_wil <- pval_comb_block(
    Z = Z, Y = Y, k = k_C, c = c_val,
    block = block,
    methods.list.all = methods.list.cmrss.wil,
    stat.null = cms_null_wil, statistic = FALSE
  )

  # Here three null draws have the same rank sum as the smallest rank sum
  # under the hypothesis, 598. RIQITE's Wilcoxon scores are the ranks
  # divided by n, which have no exact binary form, so two of the three
  # draws come out as 11.959999999999999 against a smallest value of
  # 11.960000000000001. RIQITE compares them with a plain >= and returns
  # 770/1000. CMRSS counts all three
  # (R/ties.R). Wilcoxon scores are the ranks 1 to n, so every rank sum is
  # a whole number and the p-value can be computed exactly; CMRSS is
  # checked against that count rather than against RIQITE.
  rank_sum_null <- colSums((1:n) * Z.perm)
  rank_sum_min <- min_stat(Z, Y, k_R, c_val, score = 1:n)
  expect_identical(as.numeric(cmrss_wil), sum(rank_sum_null >= rank_sum_min) / nperm)

  # --- Stephenson ---
  riq_null_ste <- RIQITE::null_dist(n, m, nperm = nperm, Z.perm = Z.perm,
                                    method.list = method.list.riqite.ste)
  scores_ste   <- list(score_all_blocks(block.sum$nb, methods.list.cmrss.ste[[1]]))
  cms_null_ste <- com_null_dist_block(Z, block, methods.list.cmrss.ste, scores_ste,
                                      null.max = nperm, weight, block.sum,
                                      Z.perm = Z.perm)

  riqite_ste <- RIQITE::pval_quantile(
    Z = Z, Y = Y, k = k_R, c = c_val,
    method.list = method.list.riqite.ste,
    stat.null = riq_null_ste, switch = FALSE
  )

  cmrss_ste <- pval_comb_block(
    Z = Z, Y = Y, k = k_C, c = c_val,
    block = block,
    methods.list.all = methods.list.cmrss.ste,
    stat.null = cms_null_ste, statistic = FALSE
  )

  expect_equal(as.numeric(riqite_ste), as.numeric(cmrss_ste), tolerance = 1e-6)
})
