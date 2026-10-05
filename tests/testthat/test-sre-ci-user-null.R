# Tests for com_block_conf_quant_larger when the caller supplies the
# permutations (Z.perm) or the null distribution (stat.null).
#
# Why this file exists:
#
# 1. The control-unit bounds come from relabeling the experiment: units that
#    were controls become "treated" (Z <- 1 - Z) and the outcome changes sign
#    (Y <- -Y). A permutation matrix describes assignments in the original
#    labeling, so in the relabeled experiment each of its columns must be
#    relabeled too (Z.perm <- 1 - Z.perm). Passing it unchanged scores the
#    wrong units in every simulated assignment. com_conf_quant_larger_cre
#    already does the relabeling; the stratified function did not, and on
#    the electric_teachers data it returned -Inf for every control bound.
#
# 2. The critical value is the entry at position floor(L * alpha) + 1 of the
#    null distribution sorted from largest to smallest, where L must be the
#    number of simulated assignments actually drawn. The function used the
#    argument null.max for L even when the caller had supplied Z.perm or
#    stat.null of a different length. With 2000 columns and the default
#    null.max = 10^4, alpha = 0.10 picked position 1001 of 2000, so the
#    procedure rejected half of the simulated assignments: a level-0.50 test
#    reported as level 0.10.
#
# 3. A single stat.null cannot serve both halves of set = "all", because the
#    relabeled experiment has nb - mb "treated" units per stratum instead of
#    mb, so its null distribution differs whenever mb != nb / 2.
#
# The designs below use 3 of 8 units treated per stratum on purpose. With
# half of each stratum treated, Z.perm and 1 - Z.perm describe the same
# design and the first defect would be invisible.

make_sre <- function(seed = 1) {
  set.seed(seed)
  s <- 3; n_per <- 8; m_per <- 3
  N <- s * n_per
  block <- factor(rep(1:s, each = n_per))
  Z <- rep(0, N)
  for (i in 1:s) Z[sample(which(block == i), m_per)] <- 1
  # A large constant effect, so every half of the experiment should produce
  # some finite lower bounds.
  Y <- rnorm(N) + 3 * Z
  methods.list.all <- list(
    lapply(1:s, function(i) list(name = "Wilcoxon", scale = FALSE)),
    lapply(1:s, function(i) list(name = "Stephenson", s = 3, scale = FALSE))
  )
  list(Z = Z, Y = Y, block = block, s = s, N = N,
       methods.list.all = methods.list.all)
}

ci_sre <- function(d, Z, Y, set, ...) {
  com_block_conf_quant_larger(Z, Y, d$block, set = set,
                              methods.list.all = d$methods.list.all,
                              opt.method = "ILP_highs", tol = 0.05,
                              alpha = 0.1, ...)
}

test_that("set = 'control' with Z.perm equals set = 'treat' in the relabeled experiment", {
  skip_if_not(solver_available("highs"), "HiGHS not available")
  d <- make_sre()
  set.seed(10)
  Zp <- assign_block(summary_block(d$Z, d$block), 1000)

  ci_control <- ci_sre(d, d$Z, d$Y, "control", Z.perm = Zp, null.max = 1000)
  # The relabeled experiment, stated directly: controls are now treated,
  # the outcome changes sign, and so does every simulated assignment.
  ci_relabel <- ci_sre(d, 1 - d$Z, -d$Y, "treat", Z.perm = 1 - Zp,
                       null.max = 1000)

  expect_equal(ci_control, ci_relabel)
  # The symptom users saw: no finite bound at all despite an effect of 3.
  expect_true(any(is.finite(ci_control)))
})

test_that("set = 'all' with Z.perm pools treat and control halves at alpha / 2", {
  skip_if_not(solver_available("highs"), "HiGHS not available")
  d <- make_sre()
  set.seed(11)
  Zp <- assign_block(summary_block(d$Z, d$block), 1000)

  ci_all <- ci_sre(d, d$Z, d$Y, "all", Z.perm = Zp, null.max = 1000)
  half <- function(Z, Y, Zperm) {
    com_block_conf_quant_larger(Z, Y, d$block, set = "treat",
                                methods.list.all = d$methods.list.all,
                                opt.method = "ILP_highs", tol = 0.05,
                                alpha = 0.05, Z.perm = Zperm, null.max = 1000)
  }
  expected <- sort(c(half(d$Z, d$Y, Zp), half(1 - d$Z, -d$Y, 1 - Zp)))

  expect_length(ci_all, d$N)
  expect_equal(ci_all, expected)
})

test_that("the critical value uses the number of columns of Z.perm, not null.max", {
  skip_if_not(solver_available("highs"), "HiGHS not available")
  d <- make_sre()
  set.seed(12)
  Zp <- assign_block(summary_block(d$Z, d$block), 2000)

  # Same 2000 assignments; only the (ignored) null.max argument differs.
  ci_matched <- ci_sre(d, d$Z, d$Y, "treat", Z.perm = Zp, null.max = 2000)
  ci_default <- ci_sre(d, d$Z, d$Y, "treat", Z.perm = Zp)
  ci_small   <- ci_sre(d, d$Z, d$Y, "treat", Z.perm = Zp, null.max = 500)

  expect_equal(ci_default, ci_matched)
  expect_equal(ci_small, ci_matched)
})

test_that("the critical value uses the length of stat.null, not null.max", {
  skip_if_not(solver_available("highs"), "HiGHS not available")
  d <- make_sre()
  block.sum <- summary_block(d$Z, d$block)
  weight <- CMRSS:::weight_scheme(block.sum, "asymp.opt")
  scores <- lapply(d$methods.list.all,
                   function(ml) CMRSS:::score_all_blocks(block.sum$nb, ml))
  set.seed(13)
  Zp <- assign_block(block.sum, 2000)
  sn <- CMRSS:::com_null_dist_block(d$Z, block.sum$block, d$methods.list.all,
                                    scores, weight = weight,
                                    block.sum = block.sum, Z.perm = Zp)
  expect_length(sn, 2000)

  ci_matched <- ci_sre(d, d$Z, d$Y, "treat", stat.null = sn, null.max = 2000)
  ci_default <- ci_sre(d, d$Z, d$Y, "treat", stat.null = sn)

  expect_equal(ci_default, ci_matched)
})

test_that("set = 'all' refuses a single stat.null for both halves", {
  skip_if_not(solver_available("highs"), "HiGHS not available")
  d <- make_sre()
  expect_error(ci_sre(d, d$Z, d$Y, "all", stat.null = rnorm(500)),
               "stat.null")
})
