# Tests for the warning about binary and heavily tied outcomes.
#
# Why this file exists:
#
# Every rank in the package is computed with rank(..., ties.method = "first"),
# so a tie between a treated and a control unit is broken by which row comes
# first in the data. Reordering the rows, with each unit keeping its own
# treatment and outcome, then changes p-values and confidence bounds. With a
# 0/1 outcome, 40 units and the same 2000 simulated assignments, one test's
# p-value ran from 0.39 to 0.93 across 20 row orders (issue #5). Changing the
# tie rule would move the numbers in the combined_stephenson_tests paper, so
# until that paper is published the package warns instead.
#
# The rule Jake chose on 2026-10-04: warn when the outcome takes exactly two
# values, or when one value is shared by more than 5 percent of the units.
# The warning comes from the four user-facing functions and from cmrss(),
# and carries the class "cmrss_ties_warning" so callers can muffle it.

set.seed(31)
n <- 40; m <- 20
Z <- sample(rep(c(1, 0), c(m, n - m)))
Y_binary <- rbinom(n, 1, 0.5)
Y_continuous <- rnorm(n)
# 3 of 40 units (7.5 percent) share the value 0; every other value is unique.
Y_tied <- c(rep(0, 3), rnorm(n - 3, mean = 5))
# 2 of 40 units (5 percent) share a value: at the threshold, not above it.
Y_two_tied <- c(rep(0, 2), rnorm(n - 2, mean = 5))
ml <- list(list(name = "Wilcoxon", scale = FALSE))

test_that("the check flags a two-valued outcome", {
  expect_warning(CMRSS:::check_outcome_ties(Y_binary),
                 class = "cmrss_ties_warning")
  expect_warning(CMRSS:::check_outcome_ties(Y_binary), "two values")
})

test_that("the check flags one value shared by more than 5 percent of units", {
  expect_warning(CMRSS:::check_outcome_ties(Y_tied),
                 class = "cmrss_ties_warning")
  # The message reports the count and the share so the user can judge.
  expect_warning(CMRSS:::check_outcome_ties(Y_tied), "3 of 40")
})

test_that("the check is silent for continuous outcomes and at the threshold", {
  expect_silent(CMRSS:::check_outcome_ties(Y_continuous))
  expect_silent(CMRSS:::check_outcome_ties(Y_two_tied))
})

test_that("a small sample with no ties is silent", {
  # In 12 units one unit is 8 percent of the sample, but a value held by one
  # unit is not a tie, and ranks cannot depend on row order.
  expect_silent(CMRSS:::check_outcome_ties(rnorm(12)))
})

test_that("the package's own example data trip the heavy-ties rule", {
  # 29 of 233 teachers share one gain score, 12 percent of the sample.
  data(electric_teachers, package = "CMRSS", envir = environment())
  expect_warning(CMRSS:::check_outcome_ties(electric_teachers$gain),
                 "29 of 233")
})

test_that("comb_p_val_cre warns on a binary outcome", {
  expect_warning(comb_p_val_cre(Z, Y_binary, k = 30, c = 0, ml, nperm = 200),
                 class = "cmrss_ties_warning")
})

test_that("com_conf_quant_larger_cre warns on a binary outcome", {
  expect_warning(com_conf_quant_larger_cre(Z, Y_binary, ml, nperm = 200),
                 class = "cmrss_ties_warning")
})

test_that("the stratified functions warn on a binary outcome", {
  skip_if_not(solver_available("highs"), "HiGHS not available")
  block <- factor(rep(1:2, each = 20))
  Zb <- unlist(lapply(1:2, function(i) sample(rep(c(1, 0), 10))))
  mla <- list(lapply(1:2, function(i) list(name = "Wilcoxon", scale = FALSE)))
  expect_warning(pval_comb_block(Zb, Y_binary, k = 15, c = 0, block, mla,
                                 opt.method = "ILP_highs", null.max = 200),
                 class = "cmrss_ties_warning")
  expect_warning(com_block_conf_quant_larger(Zb, Y_binary, block, set = "treat",
                                             methods.list.all = mla,
                                             opt.method = "ILP_highs",
                                             null.max = 200, tol = 0.1),
                 class = "cmrss_ties_warning")
})

test_that("the existing functions stay silent on a continuous outcome", {
  expect_no_warning(comb_p_val_cre(Z, Y_continuous, k = 30, c = 0, ml,
                                   nperm = 200))
})
