# CMRSS 0.2.21 gave the variables that hold the statistics (eta_h) and
# their largest standardized value (theta) finite bounds in every HiGHS
# solve of the combined stratified program, because without them HiGHS
# 1.14.0.2 had stopped with a C++ error on one program under every setting
# tried. The bounds exclude no choice of the x's, yet on other programs
# they led HiGHS to return a solution that is not the minimum with the
# status "optimal". In a simulation in 20 strata of 50 units
# (combined_stephenson_tests, code/simulations/zeta_regret/README.md,
# section 8.8) that happened in 13 of 6,120,000 p-values by more than 10
# of the 2000 null draws. The largest: replication 1690 of the cell where
# the gap between strata widens, where HiGHS returned 0.9009 with the
# bounds and 0.2546 without them, and the p-value fell from 0.5420 (the
# value Gurobi gives) to 0.2830. A minimum that is too large gives a
# p-value that is too small, so the error makes rejection easier.
#
# So, from 0.2.22, HiGHS solves the program as CMRSS 0.2.20 did, with no
# bounds, first with presolve on and then with presolve off, and the
# bounded program is solved only when both of those fail. A solve fails
# when it stops with an error or ends with a status other than optimal.
# Every program that 0.2.20 solved to an optimum gets the same answer.

## Whether a call to highs::highs_solve() carries the bounds: theta, the
## last variable, has a finite lower bound only in the bounded program.
bounded_call <- function(lower) is.finite(lower[length(lower)])

## A stand-in for highs::highs_solve() that records, for every call,
## whether it was bounded and its presolve setting, fails as `fail` says,
## and otherwise calls the real solver.
recording_highs <- function(calls, fail = function(bounded, presolve) FALSE) {
  real <- highs::highs_solve
  function(..., lower, control) {
    bounded <- bounded_call(lower)
    presolve <- if (is.null(control$presolve)) "choose" else control$presolve
    calls$log <- c(calls$log, paste(if (bounded) "bounded" else "unbounded", presolve))
    if (fail(bounded, presolve)) stop("vector::reserve")
    real(..., lower = lower, control = control)
  }
}

## The arguments pval_comb_block() (comb.method = 1) passes to
## HiGHS_sol_com() for a small stratified experiment.
small_program <- function(k = 8) {
  set.seed(20261019)
  block <- factor(rep(1:4, each = 8))
  Z <- unlist(lapply(1:4, function(b) sample(rep(0:1, each = 4))))
  Y <- rnorm(32) + Z
  ml <- CMRSS:::polynomial_methods(c(2, 6), as.vector(table(block)))
  bs <- CMRSS:::summary_block(Z, block)
  weight <- CMRSS:::weight_scheme(bs, "asymp.opt")
  scores <- lapply(ml, function(m) CMRSS:::score_all_blocks(bs$nb, m))
  list(Z = Z, block = bs$block, weight = weight,
       coeflists = CMRSS:::comb_matrix_block(Z, Y, bs$block, 0, ml, scores, bs),
       p = sum(Z) - k, ms_list = CMRSS:::mu_sigma_list(Z, bs$block, weight, ml, scores, bs),
       exact = TRUE, block.sum = bs)
}

test_that("the program on which the bounds led HiGHS to a wrong optimum solves to its minimum", {
  skip_if_not(solver_available("highs"), "HiGHS not available")
  ## The arguments of HiGHS_sol_com() in replication 1690 of the cell in 20
  ## strata of 50 where the gap between strata widens: k = 450 of the 500
  ## treated, four polynomial degrees (rule H4-f4-log), saved on koelsch.
  ## With 0.2.21's bounds HiGHS returns 0.9009211853 with status optimal;
  ## without them, and with Gurobi's p-value, the minimum is 0.2545557770.
  a <- readRDS(test_path("fixtures", "highs_wrong_optimum.rds"))
  result <- do.call(CMRSS:::HiGHS_sol_com, a)
  expect_equal(result$obj, 0.2545557770, tolerance = 1e-8)
})

test_that("a program that solves without bounds is solved once, without bounds", {
  skip_if_not(solver_available("highs"), "HiGHS not available")
  calls <- new.env(); calls$log <- character(0)
  local_mocked_bindings(highs_solve = recording_highs(calls), .package = "highs")
  do.call(CMRSS:::HiGHS_sol_com, small_program())
  expect_equal(calls$log, "unbounded choose")
})

test_that("the bounded program is solved only after both unbounded solves stop with an error, and gives the same minimum", {
  skip_if_not(solver_available("highs"), "HiGHS not available")
  a <- small_program()
  clean <- do.call(CMRSS:::HiGHS_sol_com, a)$obj
  calls <- new.env(); calls$log <- character(0)
  local_mocked_bindings(highs_solve = recording_highs(calls, fail = function(bounded, presolve) !bounded),
                        .package = "highs")
  expect_equal(do.call(CMRSS:::HiGHS_sol_com, a)$obj, clean, tolerance = 1e-9)
  expect_equal(calls$log, c("unbounded choose", "unbounded off", "bounded choose"))
})

test_that("a solve that ends without an optimum counts as a failure, and the next solve is tried", {
  skip_if_not(solver_available("highs"), "HiGHS not available")
  a <- small_program()
  clean <- do.call(CMRSS:::HiGHS_sol_com, a)$obj
  ## the unbounded solves return HiGHS's "Solve error" status with an
  ## objective of 0, as HiGHS 1.14.0.2 did on nearly fixed variables
  real <- highs::highs_solve
  calls <- new.env(); calls$log <- character(0)
  local_mocked_bindings(highs_solve = function(..., lower, control) {
    bounded <- bounded_call(lower)
    calls$log <- c(calls$log, if (bounded) "bounded" else "unbounded")
    if (!bounded) return(list(status = 4L, status_message = "Solve error", objective_value = 0,
                              primal_solution = numeric(0)))
    real(..., lower = lower, control = control)
  }, .package = "highs")
  expect_equal(do.call(CMRSS:::HiGHS_sol_com, a)$obj, clean, tolerance = 1e-9)
  expect_equal(calls$log, c("unbounded", "unbounded", "bounded"))
})

test_that("when every solve fails, CMRSS stops and reports each solve's error", {
  skip_if_not(solver_available("highs"), "HiGHS not available")
  local_mocked_bindings(highs_solve = function(..., lower, control) {
    stop(if (bounded_call(lower)) "bounded failure" else "unbounded failure")
  }, .package = "highs")
  expect_error(do.call(CMRSS:::HiGHS_sol_com, small_program()),
               "unbounded failure.*unbounded failure.*bounded failure")
})
