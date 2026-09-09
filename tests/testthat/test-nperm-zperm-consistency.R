################################################################################
# When a permutation matrix is supplied, it decides how many draws there are
#
# The functions here take both a permutation matrix Z.perm and a count nperm.
# The two say the same thing, so a caller who supplies Z.perm has already
# fixed the number of draws and whatever is in nperm is redundant. Most of the
# package already treats it that way: null_dist_multiple() sets
# nperm <- ncol(Z.perm), and com_conf_quant_larger_cre() does the same before
# it starts.
#
# comb_null_dist_cre() did not. It allocated its tail-probability matrix at
# ncol = nperm while filling it from a null distribution ncol(Z.perm) wide, so
# R recycled the values to fit. The p-value that came out of comb_p_val_cre()
# then depended on an argument that should have had no effect. On a 20-unit
# example with 200 permutations, nperm = 200 gave p = 0.8 and nperm = 1000 gave
# p = 0.000000: the same permutations, the same data, and a rejection at any
# level from an argument the caller could reasonably think was ignored.
#
# comb_p_val_cre() is exported, so this was reachable from outside the package.
# Fixed in 0.2.10 by taking the count from the null distribution actually in
# hand. These tests hold that in place.
################################################################################

## A small completely randomized experiment, fixed so the tests do not drift.
fixture <- function() {
  set.seed(1)
  n <- 20L
  m <- 12L
  Z <- rep(0, n); Z[sample(n, m)] <- 1
  list(n = n, m = m, Z = Z, Y = stats::rnorm(n),
       methods.list = list(list(name = "Wilcoxon", scale = FALSE),
                           list(name = "Stephenson", s = 3, scale = FALSE)))
}


test_that("comb_null_dist_cre returns one value per supplied permutation", {
  f <- fixture()
  set.seed(2)
  Zp <- assign_CRE(f$n, f$m, 200)

  ## The null distribution has as many entries as there are draws, and there
  ## are ncol(Zp) draws. What nperm says is beside the point.
  for (np in c(200, 1000, 5000)) {
    r <- comb_null_dist_cre(f$n, f$m, f$methods.list, Z.perm = Zp, nperm = np)
    expect_length(r, ncol(Zp))
  }
})


test_that("comb_null_dist_cre gives the same answer whatever nperm says", {
  f <- fixture()
  set.seed(2)
  Zp <- assign_CRE(f$n, f$m, 200)

  base <- comb_null_dist_cre(f$n, f$m, f$methods.list, Z.perm = Zp, nperm = ncol(Zp))
  for (np in c(1000, 5000)) {
    expect_identical(comb_null_dist_cre(f$n, f$m, f$methods.list, Z.perm = Zp, nperm = np),
                     base, info = paste("nperm =", np))
  }
})


test_that("comb_p_val_cre gives the same p-value whatever nperm says", {
  ## This is the exported function, and the one where the old behaviour showed
  ## up as a false rejection.
  f <- fixture()
  set.seed(2)
  Zp <- assign_CRE(f$n, f$m, 200)

  base <- comb_p_val_cre(f$Z, f$Y, k = 18, c = 0, f$methods.list,
                         Z.perm = Zp, nperm = ncol(Zp))
  for (np in c(1000, 5000)) {
    expect_identical(comb_p_val_cre(f$Z, f$Y, k = 18, c = 0, f$methods.list,
                                    Z.perm = Zp, nperm = np),
                     base, info = paste("nperm =", np))
  }

  ## And the p-value is a p-value: in (0, 1].
  expect_true(base > 0 && base <= 1)
})


test_that("a supplied stat.null.mult decides the count too", {
  ## The other route in: a caller who has already built the null distribution
  ## passes it directly. Its width is then what there is to work with.
  f <- fixture()
  set.seed(2)
  Zp <- assign_CRE(f$n, f$m, 200)
  snm <- null_dist_multiple(f$n, f$m, f$methods.list, Z.perm = Zp)
  expect_equal(ncol(snm), ncol(Zp))

  r <- comb_null_dist_cre(f$n, f$m, f$methods.list, Z.perm = NULL, nperm = 7777,
                          stat.null.mult = snm)
  expect_length(r, ncol(snm))
})


test_that("without a permutation matrix, nperm still sets the number of draws", {
  ## The fix must not disturb the ordinary path, where nperm is the only thing
  ## saying how many permutations to take.
  f <- fixture()
  set.seed(5)
  r <- comb_null_dist_cre(f$n, f$m, f$methods.list, nperm = 300)
  expect_length(r, 300)
})
