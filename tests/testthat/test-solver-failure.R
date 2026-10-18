# A stratified p-value rests on the smallest value of the statistic under
# the hypothesis, which CMRSS finds by solving an integer program. HiGHS can
# finish a solve without reaching an optimum: in 0.2.21 development it
# returned the status "Solve error" with objective 0 on a program whose
# optimum was -3.266. Through 0.2.20 CMRSS warned and used the objective
# anyway, so a number that is not the minimum went into the p-value, and an
# objective above the minimum makes the p-value too small.
#
# Jake decided on 2026-10-18 what CMRSS does instead. In ordinary use it
# stops with an error that says what to try. A simulation, which must keep
# going through thousands of solves, sets
# options(CMRSS.solver_failure = "na") and gets an NA objective, and so an
# NA p-value, with a warning.

## A small stratified experiment: 3 strata of 8 units, 3 treated in each.
failure_sre <- function() {
  set.seed(20261018)
  block <- factor(rep(1:3, each = 8))
  Z <- unlist(lapply(1:3, function(b) sample(rep(0:1, c(5, 3)))))
  Y <- rnorm(24) + Z
  ## polynomial scores, which both combination methods accept
  ml <- CMRSS:::polynomial_methods(c(2, 6), rep(8, 3))
  list(Z = Z, Y = Y, block = block, ml = ml)
}

## HiGHS's reply when it finishes without an optimum, as the R package
## returns it; 7 is HiGHS's status for an optimum, 4 for "Solve error".
mock_solve_error <- function() {
  local_mocked_bindings(highs_solve = function(L, ...) {
    list(status = 4L, status_message = "Solve error", objective_value = 0,
         primal_solution = rep(0, length(L)))
  }, .package = "highs", .env = parent.frame())
}

pval_at <- function(d, comb.method) {
  pval_comb_block(d$Z, d$Y, k = 6, c = 0, block = d$block, methods.list.all = d$ml,
                  null.max = 200, comb.method = comb.method, opt.method = "ILP_highs")[["p.value"]]
}

test_that("a HiGHS solve without an optimum stops with an error that says what to try", {
  skip_if_not(solver_available("highs"), "HiGHS not available")
  d <- failure_sre()
  op <- options(CMRSS.solver_failure = NULL)
  on.exit(options(op), add = TRUE)
  mock_solve_error()
  for (comb.method in 1:2) {
    err <- tryCatch(pval_at(d, comb.method), error = function(e) conditionMessage(e))
    ## HiGHS's own status, and the remedies: the relaxed program, Gurobi,
    ## and the option a simulation uses
    expect_match(err, "Solve error", fixed = TRUE)
    expect_match(err, "LP_highs", fixed = TRUE)
    expect_match(err, "ILP_gurobi", fixed = TRUE)
    expect_match(err, "CMRSS.solver_failure", fixed = TRUE)
  }
})

test_that("with options(CMRSS.solver_failure = \"na\") the p-value is NA, with a warning", {
  skip_if_not(solver_available("highs"), "HiGHS not available")
  d <- failure_sre()
  op <- options(CMRSS.solver_failure = "na")
  on.exit(options(op), add = TRUE)
  mock_solve_error()
  for (comb.method in 1:2) {
    expect_warning(p <- pval_at(d, comb.method), "Solve error", fixed = TRUE)
    expect_true(is.na(p))
  }
})

test_that("an unknown value of the option is refused", {
  skip_if_not(solver_available("highs"), "HiGHS not available")
  d <- failure_sre()
  op <- options(CMRSS.solver_failure = "skip")
  on.exit(options(op), add = TRUE)
  mock_solve_error()
  expect_error(pval_at(d, 1), "CMRSS.solver_failure")
})

test_that("a solve that reaches its optimum is unaffected by the option", {
  skip_if_not(solver_available("highs"), "HiGHS not available")
  d <- failure_sre()
  set.seed(1); p_default <- pval_at(d, 1)
  op <- options(CMRSS.solver_failure = "na")
  on.exit(options(op), add = TRUE)
  set.seed(1); p_na <- pval_at(d, 1)
  expect_identical(p_na, p_default)
  expect_false(is.na(p_default))
})
