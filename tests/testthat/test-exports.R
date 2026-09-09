################################################################################
# Tests for the package's public interface
#
# Why these exist. Until 0.2.10 the package exported the functions that answer
# a whole question -- a p-value, a set of confidence bounds -- but kept the
# building blocks those are made of internal. A user who wants to assemble a
# procedure of their own needs the blocks: a rank score vector, the minimum of
# the statistic under a bounded null, a permutation null distribution. Without
# them the only route was CMRSS:::, which is not an interface anyone should be
# asked to rely on.
#
# The combined_stephenson_tests paper repository was doing exactly that. It
# carried a 1,731-line copy of this package's code so that its analysis scripts
# could reach these six functions, and when that copy was retired the scripts
# had to bind the six out of the namespace by hand.
#
# These tests fix the public interface in place, so that a later refactor
# cannot quietly withdraw a function someone is calling.
################################################################################

## The building blocks exported in 0.2.10.
building_blocks <- c("sort_treat", "rank_score", "min_stat",
                     "null_dist", "null_dist_multiple", "comb_null_dist_cre")

## The whole-question functions exported before that.
procedures <- c("comb_p_val_cre", "com_conf_quant_larger_cre",
                "pval_comb_block", "com_block_conf_quant_larger",
                "assign_CRE", "assign_block", "summary_block",
                "method_caughey", "method_chen_li", "method_berger_boos",
                "solve_optimization", "solver_available", "parse_opt_method",
                "get_default_solver")


test_that("the building blocks are exported", {
  ex <- getNamespaceExports("CMRSS")
  for (nm in building_blocks) {
    expect_true(nm %in% ex,
                info = paste(nm, "is no longer exported; callers outside the",
                             "package would have to use CMRSS:::"))
  }
})


test_that("the functions exported before 0.2.10 are still exported", {
  ex <- getNamespaceExports("CMRSS")
  for (nm in procedures) expect_true(nm %in% ex, info = nm)
})


test_that("each exported building block is reachable through :: and is a function", {
  ## getExportedValue() fails if the name is not exported, which is the
  ## difference between an interface and an accident of the namespace.
  for (nm in building_blocks) {
    f <- getExportedValue("CMRSS", nm)
    expect_true(is.function(f), info = nm)
    expect_identical(f, getFromNamespace(nm, "CMRSS"), info = nm)
  }
})


test_that("the exported building blocks still compute what they did when internal", {
  ## Exporting must not change behaviour. These values are what the functions
  ## returned at 0.2.9, computed here from first principles rather than copied,
  ## so the test says why each answer is right.
  set.seed(11)
  n <- 20L
  m <- 12L
  Z <- rep(0, n); Z[sample(n, m)] <- 1
  Y <- round(stats::rnorm(n), 3)

  ## Wilcoxon scores are the ranks themselves.
  expect_identical(rank_score(n, list(name = "Wilcoxon", scale = FALSE)),
                   as.numeric(seq_len(n)))

  ## Stephenson scores with parameter s are choose(rank - 1, s - 1).
  expect_identical(rank_score(n, list(name = "Stephenson", s = 3, scale = FALSE)),
                   as.numeric(choose(seq_len(n) - 1, 2)))

  ## sort_treat returns the treated indices ordered by outcome, so applying it
  ## and reading off Y must give the treated outcomes in increasing order.
  ind <- sort_treat(Y, Z)
  expect_length(ind, m)
  expect_true(all(Z[ind] == 1))
  expect_false(is.unsorted(Y[ind]))

  ## With every treated unit free to have an effect above c (k = n), min_stat
  ## exempts nobody, so the statistic is the plain rank sum of the treated
  ## units computed on Y - c * Z.
  score <- rank_score(n, list(name = "Wilcoxon", scale = FALSE))
  cc <- 0.5
  adj <- Y - cc * Z
  expected <- sum(score[rank(adj, ties.method = "first")[Z == 1]])
  expect_equal(min_stat(Z, Y, k = n, c = cc, score = score), expected)

  ## A null distribution built from a supplied permutation matrix has one entry
  ## per column, and each entry is that column's treated rank sum, so its
  ## largest possible value is the sum of the m largest scores.
  set.seed(3)
  Zp <- assign_CRE(n, m, 200)
  nd <- null_dist(n, m, score = score, Z.perm = Zp)
  expect_length(nd, ncol(Zp))
  expect_true(all(nd <= sum(sort(score, decreasing = TRUE)[seq_len(m)])))
  expect_true(all(nd >= sum(sort(score)[seq_len(m)])))

  ## null_dist_multiple stacks one such row per method, and its first row must
  ## equal the single-statistic answer for that method on the same Z.perm.
  ml <- list(list(name = "Wilcoxon", scale = FALSE),
             list(name = "Stephenson", s = 3, scale = FALSE))
  ndm <- null_dist_multiple(n, m, methods.list = ml, Z.perm = Zp)
  expect_equal(dim(ndm), c(length(ml), ncol(Zp)))
  expect_equal(as.numeric(ndm[1, ]), nd)

  ## comb_null_dist_cre returns one minimum p-value per permutation, and a
  ## p-value lies in (0, 1].
  cnd <- comb_null_dist_cre(n, m, methods.list = ml, Z.perm = Zp)
  expect_length(cnd, ncol(Zp))
  expect_true(all(cnd > 0 & cnd <= 1))
})
