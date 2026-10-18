# A p-value from CMRSS is the share of null draws at or above the smallest
# value the statistic can take under the hypothesis. A null draw can equal
# that value exactly: both are sums of the same rank scores, and a rank
# statistic takes few distinct values in a small experiment. But CMRSS
# computes the two by different routes (a sum over the observed ranks, or
# an integer program, against a sum over a permuted assignment), so two
# equal values can differ in their last digits, and through 0.2.20 the
# comparison then decided the tie by rounding. When the null draw came out
# lower it was not counted, and the p-value was too small, which is the
# direction the test's promise of a type I error rate of at most alpha in
# any sample size cannot allow. The completely randomized combined test had
# a second form of the same fault: it compared a share computed as
# 1 - (rank - 1) / nperm with one computed as a mean, and 1 - 997/1000 is
# not 3/1000 in floating point.
#
# From 0.2.21 CMRSS counts two values as equal when they differ by less
# than 1e-9 times the larger of 1 and their size. Rounding moves these
# values by about 1e-15 of their size. In each experiment below, two
# values that truly differ are at least 1e-5 apart.
#
# The tests compare CMRSS with p-values computed in whole numbers. The
# polynomial score of degree d, (r / (n + 1))^(d - 1), is a positive
# multiple of r^(d - 1), so the whole-number scores r^(d - 1) order every
# statistic the same way. With them every statistic is an exact whole
# number and every tie is decided exactly, so the p-value computed from
# them is the one CMRSS should return.

## the polynomial scores the paper uses, and whole-number scores in the
## same order
poly_method <- function(d) list(name = "Polynomial", r = d, std = TRUE, scale = FALSE)
whole_scores <- function(n, d) (1:n)^(d - 1)

## ---- completely randomized experiments

## 40 units, 20 treated, 200 null assignments
cre_experiment <- function(seed) {
  set.seed(seed)
  n <- 40; m <- 20
  Z <- sample(rep(0:1, c(n - m, m)))
  list(n = n, m = m, Z = Z, Y = rnorm(n) + Z * rexp(n),
       Z.perm = replicate(200, sample(Z)))
}

## For each degree (columns) and null draw (rows), how many null draws are
## at or above that draw, in whole numbers.
cre_null_counts <- function(e, degrees) {
  counts <- lapply(degrees, function(d) {
    x <- null_dist(e$n, e$m, score = whole_scores(e$n, d), Z.perm = e$Z.perm)
    vapply(x, function(v) sum(x >= v), numeric(1))
  })
  matrix(unlist(counts), ncol = length(degrees))
}

## For each degree, how many null draws are at or above the smallest value
## of the statistic under the hypothesis that at most n - k units have
## effects above c.
cre_observed_counts <- function(e, degrees, k, c) {
  vapply(degrees, function(d) {
    sc <- whole_scores(e$n, d)
    sum(null_dist(e$n, e$m, score = sc, Z.perm = e$Z.perm) >=
          min_stat(e$Z, e$Y, k, c, score = sc))
  }, numeric(1))
}

## The combined p-value in whole numbers: the share of null draws whose
## smallest count over the degrees is at or below the observed smallest
## count.
cre_whole_p <- function(e, degrees, k, c, null_counts = cre_null_counts(e, degrees)) {
  sum(apply(null_counts, 1, min) <= min(cre_observed_counts(e, degrees, k, c))) / ncol(e$Z.perm)
}

## ---- stratified experiments

## 4 strata of 10 units, 5 treated in each, 200 null assignments, each
## permuting the assignment within strata. With strata of equal size the
## stratified statistic is the sum of the strata's statistics, rescaled.
sre_experiment <- function(seed) {
  set.seed(seed)
  B <- 4; nb <- 10
  block <- factor(rep(1:B, each = nb))
  Z <- unlist(lapply(1:B, function(b) sample(rep(0:1, each = nb / 2))))
  Y <- rnorm(B * nb) + Z * rexp(B * nb)
  Z.perm <- replicate(200, unlist(lapply(1:B, function(b) sample(Z[block == b]))))
  list(B = B, nb = nb, mb = nb / 2, n = B * nb, m = sum(Z), block = block,
       Z = Z, Y = Y, Z.perm = Z.perm)
}

## CMRSS's null draws for the experiment: each degree's stratified
## statistic standardized, and the largest over the degrees.
sre_null <- function(e, degrees) {
  bs <- summary_block(e$Z, e$block)
  CMRSS:::com_null_dist_block(e$Z, e$block, CMRSS:::polynomial_methods(degrees, rep(e$nb, e$B)),
                              null.max = ncol(e$Z.perm), weight = CMRSS:::weight_scheme(bs),
                              block.sum = bs, Z.perm = e$Z.perm)
}

## The p-value for one degree in whole numbers, for the hypothesis that at
## most m - k treated units have effects above c. The smallest stratified
## statistic is the smallest, over the ways of spreading the m - k exempt
## units across strata, of the sum of each stratum's smallest statistic.
## The null draws take a unit's rank in its stratum to be its position
## there, as CMRSS's do.
sre_whole_p <- function(e, d, k, c) {
  sc <- whole_scores(e$nb, d)
  stratum_min <- lapply(levels(e$block), function(b) {
    i <- e$block == b
    vapply(0:e$mb, function(ii) min_stat(e$Z[i], e$Y[i], k = e$nb - ii, c = c, score = sc), numeric(1))
  })
  spread <- as.matrix(expand.grid(rep(list(0:e$mb), e$B)))
  spread <- spread[rowSums(spread) <= e$m - k, , drop = FALSE]
  S_min <- min(apply(spread, 1, function(a) sum(vapply(1:e$B, function(b) stratum_min[[b]][a[b] + 1], numeric(1)))))
  rank_in_stratum <- ave(seq_len(e$n), e$block, FUN = rank)
  S_null <- colSums(sc[rank_in_stratum] * e$Z.perm)
  sum(S_null >= S_min) / ncol(e$Z.perm)
}

## ---- the tests

test_that("values that differ only by rounding count as tied, and values that truly differ do not", {
  ## 0.1 + 0.2 is 0.30000000000000004 in floating point
  draws <- c(0.3, 1, 2)
  expect_equal(mean(draws >= 0.1 + 0.2), 2 / 3)   # the comparison through 0.2.20
  expect_equal(CMRSS:::share_at_least(draws, 0.1 + 0.2), 1)
  expect_equal(CMRSS:::share_at_most(0.1 + 0.2, 0.3), 1)
  ## a difference of a millionth is a real difference
  expect_equal(CMRSS:::share_at_least(c(0.3 - 1e-6, 1, 2), 0.3), 2 / 3)
  expect_equal(CMRSS:::share_at_most(0.3 + 1e-6, 0.3), 0)
  ## the margin grows with the size of the values
  expect_equal(CMRSS:::share_at_least(1000 - 1e-10, 1000), 1)
  expect_equal(CMRSS:::share_at_least(1000 - 1e-3, 1000), 0)
})

test_that("a single statistic's p-value in a completely randomized experiment is the one computed in whole numbers", {
  for (seed in 1:10) {
    e <- cre_experiment(seed)
    for (k in c(21, 30, 35, 40)) for (d in c(2, 4)) {
      p <- CMRSS:::pval_cre(e$Z, e$Y, k, 0, method.list = poly_method(d), Z.perm = e$Z.perm)
      expect_identical(p, cre_observed_counts(e, d, k, 0) / ncol(e$Z.perm),
                       label = sprintf("seed %d, k %d, degree %d", seed, k, d))
    }
  }
})

test_that("the combined p-value in a completely randomized experiment is the one computed in whole numbers", {
  degrees <- c(2, 4)
  for (seed in 1:10) {
    e <- cre_experiment(seed)
    null_counts <- cre_null_counts(e, degrees)
    for (k in c(21, 30, 35, 40)) {
      p <- comb_p_val_cre(e$Z, e$Y, k, 0, lapply(degrees, poly_method), Z.perm = e$Z.perm, nperm = 200)
      expect_identical(p, cre_whole_p(e, degrees, k, 0, null_counts),
                       label = sprintf("seed %d, k %d", seed, k))
    }
  }
})

test_that("a stratified p-value is the one computed in whole numbers", {
  skip_if_not(solver_available("highs"), "HiGHS not available")
  for (seed in 1:10) {
    e <- sre_experiment(seed)
    stat.null <- sre_null(e, 2)
    for (k in c(1, 8, 15, 20)) {
      p <- pval_comb_block(e$Z, e$Y, k, 0, e$block, CMRSS:::polynomial_methods(2, rep(e$nb, e$B)),
                           stat.null = stat.null, opt.method = "ILP_highs")[["p.value"]]
      expect_identical(p, sre_whole_p(e, 2, k, 0), label = sprintf("seed %d, k %d", seed, k))
    }
  }
})

test_that("when every null draw is the observed assignment, the p-value is 1", {
  ## Under the hypothesis that no treated unit's effect exceeds 0 (k = m in
  ## a stratified experiment, k = n in a completely randomized one), the
  ## smallest value of the statistic is its observed value. A null draw
  ## that repeats the observed assignment has that value too, provided a
  ## unit's position is its rank, as in CMRSS's null draws. So the outcomes
  ## are sorted (within strata), and every null draw is the observed
  ## assignment.
  skip_if_not(solver_available("highs"), "HiGHS not available")
  degrees <- c(2, 6)
  for (seed in 1:6) {
    e <- sre_experiment(seed)
    o <- order(e$block, e$Y)
    e$Z <- e$Z[o]; e$Y <- e$Y[o]; e$block <- e$block[o]
    e$Z.perm <- matrix(e$Z, e$n, 200)
    p <- pval_comb_block(e$Z, e$Y, e$m, 0, e$block, CMRSS:::polynomial_methods(degrees, rep(e$nb, e$B)),
                         stat.null = sre_null(e, degrees), opt.method = "ILP_highs")[["p.value"]]
    expect_identical(p, 1, label = sprintf("stratified, seed %d", seed))

    e <- cre_experiment(seed)
    o <- order(e$Y)
    e$Z <- e$Z[o]; e$Y <- e$Y[o]
    e$Z.perm <- matrix(e$Z, e$n, 200)
    p <- comb_p_val_cre(e$Z, e$Y, e$n, 0, lapply(degrees, poly_method), Z.perm = e$Z.perm, nperm = 200)
    expect_identical(p, 1, label = sprintf("completely randomized, seed %d", seed))
  }
})

## The interval functions find, for each k, the smallest c on a grid of
## step tol that the combined test does not reject. So at each limit L the
## whole-number p-value must exceed alpha, and at L - tol it must be at
## most alpha. Two limits need no check: -Inf, where no c was rejected, and
## a limit equal to the next k's limit, which caps it.
expect_limits_match <- function(L, k_of, whole_p, alpha, tol, what) {
  for (j in seq_along(L)) {
    if (!is.finite(L[j])) next
    capped <- j < length(L) && L[j] == L[j + 1]
    if (!capped) {
      expect_gt(whole_p(k_of(j), L[j]), alpha, label = sprintf("%s, k %d: p at its limit", what, k_of(j)))
    }
    expect_lte(whole_p(k_of(j), L[j] - tol), alpha,
               label = sprintf("%s, k %d: p just below its limit", what, k_of(j)))
  }
}

test_that("the completely randomized interval's limits agree with the p-values computed in whole numbers", {
  degrees <- c(2, 4); alpha <- 0.1; tol <- 1e-3
  for (seed in c(6, 10, 15)) {
    e <- cre_experiment(seed)
    null_counts <- cre_null_counts(e, degrees)
    L <- com_conf_quant_larger_cre(e$Z, e$Y, lapply(degrees, poly_method), nperm = 200, set = "treat",
                                   Z.perm = e$Z.perm, alpha = alpha, tol = tol)
    expect_limits_match(L, function(j) e$n - e$m + j,
                        function(k, c) cre_whole_p(e, degrees, k, c, null_counts),
                        alpha, tol, sprintf("seed %d", seed))
  }
})

test_that("the stratified interval's limits agree with the p-values computed in whole numbers", {
  skip_if_not(solver_available("highs"), "HiGHS not available")
  alpha <- 0.1; tol <- 0.01
  for (seed in 1:3) {
    e <- sre_experiment(seed)
    L <- com_block_conf_quant_larger(e$Z, e$Y, e$block, set = "treat",
                                     methods.list.all = CMRSS:::polynomial_methods(2, rep(e$nb, e$B)),
                                     opt.method = "ILP_highs", stat.null = sre_null(e, 2),
                                     tol = tol, alpha = alpha)
    ## L[j] is the limit for the j-th smallest treated effect, k = j
    expect_limits_match(L, function(j) j, function(k, c) sre_whole_p(e, 2, k, c),
                        alpha, tol, sprintf("seed %d", seed))
  }
})
