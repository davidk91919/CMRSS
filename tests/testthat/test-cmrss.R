# Tests for cmrss(), the formula and data frame interface.
#
# Why this file exists:
#
# cmrss() adds no statistics of its own. It reads the design from a formula,
# builds the nested score lists, translates a proportion into each existing
# function's k, and calls the existing functions. So the tests check that it
# returns exactly what a careful user of the existing functions would get by
# hand, with the same random seed, and that each translation is right.
#
# Decisions Jake made on 2026-10-04 that these tests pin:
#   1. Strata go after a bar in the formula: gain ~ TxAny | Site.
#   2. The user gives a proportion of a group (quantile = 0.9, set = "treat");
#      cmrss() converts it to k and reports the k it used.
#   3. The user gives Stephenson parameters as a vector s; any s larger than
#      a block's size becomes that block's size, so a small block scores 1
#      when its highest outcome is treated instead of dropping out.
#   3a. The default s is 2, round(sqrt(2 * s_max)) and s_max, where
#      s_max = min(floor(4 m / q_min), floor(n / 2)) and q_min is the smallest
#      q with (m/n)((m-1)/(n-1))...((m-q+1)/(n-q+1)) <= alpha. See the
#      section "Choosing s" in ?cmrss for the reasoning.
#   4. The tie warning (see test-outcome-ties-warning.R) comes once per call.

# A small completely randomized experiment and a small stratified one. The
# stratified design has a block of 5 units so that s = 30 must be capped.
make_cre <- function() {
  set.seed(101)
  n <- 30
  d <- data.frame(z = sample(rep(c(1, 0), c(15, 15))))
  d$y <- rnorm(n) + 2 * d$z
  d
}
make_sre <- function() {
  set.seed(102)
  sizes <- c(5, 12, 13)
  site <- factor(rep(c("a", "b", "c"), sizes))
  z <- unlist(lapply(sizes, function(nb) sample(rep(c(1, 0), c(ceiling(nb / 2), floor(nb / 2))))))
  data.frame(site = site, z = z, y = rnorm(sum(sizes)) + 2 * z)
}
# make_cre(): n = 30, m = 15 -> q_min = 4, s_max = 15, default 2, 5, 15.
s_default <- c(2, 5, 15)
# make_sre(): n = 30, m = 16 -> q_min = 5 (probability 0.0307),
# s_max = min(floor(64 / 5), 15) = 12, default 2, 5, 12.
s_sre_default <- c(2, 5, 12)
steph <- function(s) lapply(s, function(x) list(name = "Stephenson", s = x, scale = FALSE))
# The default scores are polynomial, (r / (n + 1))^(s - 1); see ?cmrss.
poly <- function(s) lapply(s, function(x) list(name = "Polynomial", r = x, std = TRUE, scale = FALSE))

## The default s ----------------------------------------------------------------

test_that("q_min is the fewest treated units that could reach alpha", {
  # n = 200, m = 100: 4 units give probability 0.0606, 5 give 0.0297.
  expect_equal(CMRSS:::min_detectable_treated(200, 100, alpha = 0.05), 5)
  expect_equal(prod((100 - 0:3) / (200 - 0:3)), 0.0606, tolerance = 1e-3)
  expect_equal(prod((100 - 0:4) / (200 - 0:4)), 0.0297, tolerance = 1e-3)
})

test_that("the default s matches the table in ?cmrss", {
  ds <- function(n, m) CMRSS:::default_score_parameters(n, m, alpha = 0.05)
  expect_equal(ds(200, 100), c(2, 13, 80))
  expect_equal(ds(233, 164), c(2, 12, 72))  # electric_teachers
  expect_equal(ds(30, 15), c(2, 5, 15))
  expect_equal(ds(100, 20), c(2, 9, 40))
})

test_that("the default s never exceeds n / 2", {
  # n = 1000, m = 500: q_min = 5 and 4 m / q_min = 400 is below n / 2 = 500,
  # so 400 stays. n = 40, m = 36: whatever q_min is, s_max must be <= 20.
  expect_lte(max(CMRSS:::default_score_parameters(40, 36, alpha = 0.05)), 20)
  expect_equal(max(CMRSS:::default_score_parameters(1000, 500, alpha = 0.05)), 400)
})

test_that("a design too small to reach alpha is refused", {
  # n = 4, m = 2: even both treated units holding the top two ranks has
  # probability (2/4)(1/3) = 0.167 > 0.05, and q cannot exceed m.
  expect_error(CMRSS:::default_score_parameters(4, 2, alpha = 0.05), "alpha")
})

## Building the score lists ---------------------------------------------------

test_that("an s larger than a block's size becomes that block's size", {
  ml <- CMRSS:::stephenson_methods(s = c(2, 30), nb = c(5, 40))
  # One element per statistic, each with one specification per block.
  expect_length(ml, 2)
  expect_length(ml[[2]], 2)
  expect_equal(ml[[2]][[1]]$s, 5)   # capped in the block of 5
  expect_equal(ml[[2]][[2]]$s, 30)  # unchanged in the block of 40
  expect_equal(ml[[1]][[1]]$s, 2)
})

test_that("a capped block scores 1 for its top rank and 0 elsewhere", {
  # choose(r - 1, 4) for r = 1..5 is 0, 0, 0, 0, 1: the score marks the
  # block's highest outcome, so the block does not drop out.
  ml <- CMRSS:::stephenson_methods(s = 30, nb = 5)
  expect_equal(rank_score(5, ml[[1]][[1]]), c(0, 0, 0, 0, 1))
})

test_that("polynomial parameters are not lowered in small blocks", {
  # A polynomial score is never 0, so a block of 5 keeps zeta = 30; the
  # package's Polynomial score with std = TRUE is (r / (n_b + 1))^(zeta - 1).
  ml <- CMRSS:::polynomial_methods(zeta = c(2, 30), nb = c(5, 40))
  expect_equal(ml[[2]][[1]]$r, 30)
  expect_true(isTRUE(ml[[2]][[1]]$std))
  expect_equal(rank_score(5, ml[[2]][[1]]), ((1:5) / 6)^29)
})

## Completely randomized experiments ------------------------------------------

test_that("CRE bounds equal com_conf_quant_larger_cre with the same seed", {
  d <- make_cre()
  set.seed(1)
  fit <- cmrss(y ~ z, data = d, set = "treat", nperm = 500, alpha = 0.1)
  set.seed(1)
  by_hand <- com_conf_quant_larger_cre(d$z, d$y, poly(s_default), nperm = 500,
                                       set = "treat", alpha = 0.1, tol = 0.01)
  expect_equal(fit$bounds$lower, by_hand)
  expect_equal(fit$bounds$k, 1:15)
  expect_equal(fit$bounds$proportion, (1:15) / 15)
})

test_that("CRE test: quantile 0.9 of the treated becomes k = n - m + floor(0.9 m)", {
  d <- make_cre()
  set.seed(2)
  fit <- cmrss(y ~ z, data = d, quantile = 0.9, c = 1, set = "treat",
               nperm = 500)
  set.seed(2)
  # comb_p_val_cre counts k over all n units: the floor(0.9 * 15) = 13th
  # smallest treated effect is the (30 - 15 + 13) = 28th in its numbering.
  p_by_hand <- comb_p_val_cre(d$z, d$y, k = 28, c = 1, poly(s_default),
                              nperm = 500)
  expect_equal(fit$test$k, 13)
  expect_equal(fit$test$p.value, p_by_hand)
})

test_that("a minus sign on the outcome tests effects on -y", {
  # The harm question uses effects on -y, which are minus the effects on y.
  d <- make_cre()
  set.seed(3)
  fit <- cmrss(-y ~ z, data = d, set = "treat", nperm = 500)
  set.seed(3)
  by_hand <- com_conf_quant_larger_cre(d$z, -d$y, poly(s_default),
                                       nperm = 500, set = "treat",
                                       alpha = 0.05, tol = 0.01)
  expect_equal(fit$bounds$lower, by_hand)
})

## Stratified experiments -------------------------------------------------------

test_that("SRE bounds equal com_block_conf_quant_larger with the default polynomial scores", {
  skip_if_not(solver_available("highs"), "HiGHS not available")
  d <- make_sre()
  set.seed(4)
  fit <- cmrss(y ~ z | site, data = d, set = "treat", nperm = 300,
               tol = 0.05, opt.method = "ILP_highs")
  set.seed(4)
  ml <- CMRSS:::polynomial_methods(s_sre_default, nb = c(5, 12, 13))
  by_hand <- com_block_conf_quant_larger(d$z, d$y, d$site, set = "treat",
                                         methods.list.all = ml,
                                         opt.method = "ILP_highs",
                                         null.max = 300, tol = 0.05,
                                         alpha = 0.05)
  expect_equal(fit$bounds$lower, by_hand)
})

test_that("SRE test: quantile 0.9 of the treated becomes k = floor(0.9 m)", {
  skip_if_not(solver_available("highs"), "HiGHS not available")
  d <- make_sre()
  m <- sum(d$z)
  set.seed(5)
  fit <- cmrss(y ~ z | site, data = d, quantile = 0.9, c = 1, set = "treat",
               nperm = 300, opt.method = "ILP_highs")
  set.seed(5)
  ml <- CMRSS:::polynomial_methods(s_sre_default, nb = c(5, 12, 13))
  # pval_comb_block counts k over the treated units only.
  p_by_hand <- pval_comb_block(d$z, d$y, k = floor(0.9 * m), c = 1, d$site,
                               ml, null.max = 300, opt.method = "ILP_highs",
                               statistic = FALSE)
  expect_equal(fit$test$k, floor(0.9 * m))
  expect_equal(fit$test$p.value, p_by_hand)
})

test_that("scores = 'stephenson' uses the default values as s, lowered in small blocks", {
  skip_if_not(solver_available("highs"), "HiGHS not available")
  d <- make_sre()
  set.seed(7)
  fit <- cmrss(y ~ z | site, data = d, set = "treat", scores = "stephenson",
               nperm = 300, tol = 0.05, opt.method = "ILP_highs")
  set.seed(7)
  ml <- CMRSS:::stephenson_methods(s_sre_default, nb = c(5, 12, 13))
  by_hand <- com_block_conf_quant_larger(d$z, d$y, d$site, set = "treat",
                                         methods.list.all = ml,
                                         opt.method = "ILP_highs",
                                         null.max = 300, tol = 0.05,
                                         alpha = 0.05)
  expect_equal(fit$bounds$lower, by_hand)
})

## Inputs ------------------------------------------------------------------------

test_that("a treatment that is not 0/1 or logical is refused with a clear message", {
  d <- make_cre()
  d$z2 <- d$z + 1
  expect_error(cmrss(y ~ z2, data = d, nperm = 100), "0 and 1")
})

test_that("missing values are refused and counted", {
  d <- make_cre()
  d$y[c(2, 5)] <- NA
  expect_error(cmrss(y ~ z, data = d, nperm = 100), "2 rows")
})

test_that("quantile must be a proportion strictly between 0 and 1", {
  d <- make_cre()
  expect_error(cmrss(y ~ z, data = d, quantile = 1.5, nperm = 100), "quantile")
})

## Output ------------------------------------------------------------------------

test_that("print reports how many units have effects above c", {
  d <- make_cre()
  set.seed(6)
  fit <- cmrss(y ~ z, data = d, set = "treat", nperm = 500)
  n_above <- sum(fit$bounds$lower > 0)
  expect_output(print(fit), paste0("at least ", n_above, " of 15"))
})

test_that("print gives the confidence level without rounding", {
  # alpha = 0.025, used when two analyses are reported together, is 97.5
  # percent confidence, not 98.
  d <- make_cre()
  set.seed(9)
  fit <- cmrss(y ~ z, data = d, set = "treat", alpha = 0.025, nperm = 200,
               tol = 0.1)
  expect_output(print(fit), "97.5 percent confidence")
})

test_that("print says which outcome the effects are on", {
  # With -y on the left the count is of effects on -y, that is, units whose
  # effect on y is negative; the printed line must not hide the sign.
  d <- make_cre()
  set.seed(8)
  fit <- cmrss(-y ~ z, data = d, set = "treat", nperm = 200, tol = 0.1)
  expect_output(print(fit), "effects on -y above 0")
})

test_that("a tied outcome produces one warning per call, not several", {
  d <- make_cre()
  d$y <- round(d$y)  # few distinct values
  w <- 0
  withCallingHandlers(
    cmrss(y ~ z, data = d, set = "treat", nperm = 200),
    cmrss_ties_warning = function(cnd) {
      w <<- w + 1
      invokeRestart("muffleWarning")
    }
  )
  expect_equal(w, 1)
})
