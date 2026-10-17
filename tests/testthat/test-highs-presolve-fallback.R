# HiGHS's presolve, the pass in which it simplifies an integer program
# before solving it, sometimes stops with the C++ error "vector::reserve".
# On the Illinois Campus Cluster (Linux, x86-64, highs 1.14.0.2) it did so
# in 16 of 180,000 replications of a simulation in 20 strata of 50 units
# (combined_stephenson_tests, code/simulations/zeta_regret/README.md,
# section 8.8). Each time the stratified test stopped, so a researcher
# without Gurobi, whom CMRSS sends to HiGHS, could not get a p-value. The
# saved program that failed solves on the cluster with presolve off, and on
# the laptop with presolve on, to the same optimum, -5.705271041. So when a
# HiGHS solve stops with an error, CMRSS solves the same program again with
# presolve off. The second solve runs only after an error, so no p-value
# that CMRSS computed before this change can change.
#
# The error depends on the machine and on the last digits of the program,
# so the first two tests make HiGHS fail on purpose: a stand-in for
# highs::highs_solve() throws the cluster's error whenever presolve is on
# and otherwise calls the real solver.

## A stratified experiment small enough to solve in a second, with two
## statistics combined, as in the paper's stratified test.
small_sre <- function() {
  set.seed(20261017)
  block <- factor(rep(1:4, each = 10))
  Z <- unlist(lapply(1:4, function(b) sample(rep(0:1, each = 5))))
  Y <- rnorm(40) + Z
  nb <- as.vector(table(block))
  list(Z = Z, Y = Y, block = block, ml = CMRSS:::polynomial_methods(c(2, 6), nb))
}

## The stand-in solver: it fails like the cluster's HiGHS while presolve
## is on, records every call's presolve setting, and otherwise solves.
failing_highs <- function(calls) {
  real <- highs::highs_solve
  function(..., control) {
    calls$presolve <- c(calls$presolve, if (is.null(control$presolve)) "choose" else control$presolve)
    if (!identical(control$presolve, "off")) stop("vector::reserve")
    real(..., control = control)
  }
}

test_that("a HiGHS solve that stops with an error is solved again with presolve off, and the p-value is unchanged", {
  skip_if_not(solver_available("highs"), "HiGHS not available")
  d <- small_sre()
  ## both combination methods: comb.method 1 solves the combined program
  ## (HiGHS_sol_com), 2 the stratum-by-stratum one (HiGHS_sol_stratum_com)
  pval <- function(comb) {
    set.seed(1)
    pval_comb_block(d$Z, d$Y, k = 15, c = 0, block = d$block, methods.list.all = d$ml,
                    null.max = 200, opt.method = "ILP_highs", comb.method = comb)
  }
  p_clean <- c(pval(1), pval(2))
  calls <- new.env(); calls$presolve <- character(0)
  local_mocked_bindings(highs_solve = failing_highs(calls), .package = "highs")
  p_fallback <- c(pval(1), pval(2))
  expect_equal(p_fallback, p_clean)
  ## every solve failed once with presolve on and was solved again with it off
  expect_gt(length(calls$presolve), 0)
  expect_equal(calls$presolve, rep(c("choose", "off"), length(calls$presolve) / 2))
})

test_that("an error that persists with presolve off stops the test and reports both errors", {
  skip_if_not(solver_available("highs"), "HiGHS not available")
  d <- small_sre()
  local_mocked_bindings(highs_solve = function(..., control) {
    stop(if (identical(control$presolve, "off")) "second failure" else "vector::reserve")
  }, .package = "highs")
  expect_error(pval_comb_block(d$Z, d$Y, k = 15, c = 0, block = d$block, methods.list.all = d$ml,
                               null.max = 200, opt.method = "ILP_highs"),
               "vector::reserve.*second failure")
})

test_that("the program that stopped HiGHS on the Campus Cluster solves to its optimum", {
  skip_if_not(solver_available("highs"), "HiGHS not available")
  ## 523 variables and 25 constraints: the combined test's program for
  ## replication 1302 of one cell, saved on the cluster just before HiGHS
  ## failed on it. On the cluster's build the first solve fails and the
  ## second succeeds; on other machines the first solve succeeds.
  args <- readRDS(test_path("fixtures", "highs_presolve_failure.rds"))
  result <- do.call(CMRSS:::highs_solve_presolve_fallback, args)
  expect_equal(result$status, 7)   # 7 = optimal in HiGHS
  expect_equal(result$objective_value, -5.705271041, tolerance = 1e-9)
})
