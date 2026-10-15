# Tests for divide_by_top, the Stephenson score divided by the top unit's
# score, written before the code.
#
# Why: rank_score() computes the Stephenson score choose(r - 1, s - 1)
# directly. Once n reaches 1,032, choose(n - 1, 499) exceeds the largest
# double, so with s = 500 the top scores are Inf and every p-value is NA.
# The regret study in combined_stephenson_tests hit this on its pattern
# with 1,157 units.
#
# The fix, decided with Jake on 2026-10-14: an argument divide_by_top in the
# Stephenson method list. TRUE, and also when the argument is left out,
# divides each score by the top unit's score, choose(n - 1, s - 1),
# computed on the log scale as exp(lchoose(r - 1, s - 1) - lchoose(n - 1,
# s - 1)), which cannot overflow. FALSE gives the original binomial
# coefficients. In an experiment without strata the two give the same
# p-values: the observed statistic and every null draw are multiplied by
# the same constant, so the share of null draws at or above the observed
# statistic does not change. With strata of unequal size they differ,
# because each stratum's scores are divided by a different number.

stephenson <- function(s, ...) list(name = "Stephenson", s = s, scale = FALSE, ...)

test_that("divide_by_top = TRUE gives the raw scores divided by the top score", {
  n <- 30; s <- 5
  raw <- choose(0:(n - 1), s - 1)
  divided <- rank_score(n, stephenson(s, divide_by_top = TRUE))
  expect_equal(divided, raw / raw[n])
  # ranks below s score exactly 0 and the top rank scores exactly 1
  expect_true(all(divided[1:(s - 1)] == 0))
  expect_identical(divided[n], 1)
})

test_that("leaving divide_by_top out divides, and FALSE keeps the binomial coefficients", {
  n <- 30; s <- 5
  expect_identical(rank_score(n, stephenson(s)),
                   rank_score(n, stephenson(s, divide_by_top = TRUE)))
  expect_identical(rank_score(n, stephenson(s, divide_by_top = FALSE)),
                   choose(0:(n - 1), s - 1))
})

test_that("divided Stephenson scores are finite where the raw ones overflow", {
  n <- 1157; s <- 500
  # the failure this fixes: the raw top score is Inf
  expect_true(is.infinite(choose(n - 1, s - 1)))
  divided <- rank_score(n, stephenson(s))
  expect_true(all(is.finite(divided)))
  expect_true(all(divided >= 0 & divided <= 1))
  expect_false(is.unsorted(divided))
  expect_equal(divided[n], 1)
  # choose(r - 1, s - 1) / choose(r - 2, s - 1) = (r - 1) / (r - s), so the
  # second-highest score is (n - s) / (n - 1) of the top one
  expect_equal(divided[n - 1] / divided[n], (n - s) / (n - 1), tolerance = 1e-12)
})

test_that("without strata, divided and raw Stephenson scores give the same p-values", {
  set.seed(20261014)
  n <- 60; m <- 30
  Z <- sample(rep(c(1, 0), c(m, n - m)))
  Y <- rnorm(n) + Z * rexp(n)
  Z.perm <- assign_CRE(n, m, 2000)
  raw_methods <- lapply(c(2, 6, 20), function(s) stephenson(s, divide_by_top = FALSE))
  divided_methods <- lapply(c(2, 6, 20), function(s) stephenson(s, divide_by_top = TRUE))
  # k = n is the null hypothesis of no positive effect on any unit; the
  # smaller k leave n - k units free to take any effect
  for (k in c(n, n - 5, n - 15)) {
    expect_equal(comb_p_val_cre(Z, Y, k = k, c = 0, methods.list = divided_methods, Z.perm = Z.perm, nperm = 2000),
                 comb_p_val_cre(Z, Y, k = k, c = 0, methods.list = raw_methods, Z.perm = Z.perm, nperm = 2000))
  }
})
