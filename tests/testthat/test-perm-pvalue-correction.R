# Tests for the Monte Carlo permutation p-value correction of
# Phipson and Smyth (2010), "Permutation P-values Should Never Be Zero:
# Calculating Exact P-values When Permutations Are Randomly Drawn,"
# Statistical Applications in Genetics and Molecular Biology 9(1),
# doi:10.2202/1544-6115.1585.  A copy is at
# references/Phipson_Smyth_2010_permutation_pvalues.pdf.
#
# The problem in CMRSS.  `pval_comb_block` computes
#
#     pval <- mean(stat.null >= stat.min)
#
# at R/CMRSS_SRE.R:54 and :71.  `stat.null` holds `null.max` draws from the
# reference distribution, so this is the proportion of draws at least as
# extreme as the observed statistic.  When no draw reaches `stat.min` the
# result is exactly 0, which claims a p-value smaller than the simulation
# can resolve.  Phipson and Smyth show that substituting such an unbiased
# estimator for the exact p-value inflates the type I error rate, and that
# adding one to both numerator and denominator does not.
#
# Why the corrected form is valid here.  Under H_{k,c}^treat the quantity
# tilde_t(Z, Y - Z*xi) has exactly the reference distribution, because the
# reference distribution is distribution-free (main.tex Theorem at line 900:
# G depends on the design and the scores, not on the outcomes).  The
# statistic CMRSS actually compares against it, `stat.min`, is a minimum
# over the null set and so is no larger than that quantity.  The observed
# value is therefore stochastically no larger than a draw, and counting it
# among the draws gives a p-value that is valid and, if anything,
# conservative.
#
# These tests describe a helper that does not exist yet.  They skip until
# it does.  See PLAN.md item 4A: wiring the correction into
# `pval_comb_block` changes every p-value by roughly 1/null.max, which is a
# numerical-output change and waits on the R&R replication run.

has_helper <- function() {
  exists("perm_pvalue", where = asNamespace("CMRSS"), inherits = FALSE)
}

skip_unless_helper <- function() {
  skip_if_not(has_helper(),
              "perm_pvalue() not implemented yet (PLAN.md item 4A)")
}


test_that("perm_pvalue matches the (b + 1) / (m + 1) formula", {
  skip_unless_helper()

  # A reference distribution we can count by hand: the integers 1..100.
  stat.null <- as.numeric(1:100)

  # 10 of the 100 draws are >= 91, so b = 10 and m = 100.
  expect_equal(CMRSS:::perm_pvalue(stat.null, 91), 11 / 101)

  # Every draw is >= 1, so b = m = 100 and the p-value is 1.
  expect_equal(CMRSS:::perm_pvalue(stat.null, 1), 101 / 101)

  # No draw reaches 101, so b = 0 and the p-value is the smallest the
  # simulation can report rather than 0.
  expect_equal(CMRSS:::perm_pvalue(stat.null, 101), 1 / 101)
})


test_that("perm_pvalue never returns zero", {
  skip_unless_helper()

  # The uncorrected form returns 0 here; that is the failure the paper is
  # named for.
  stat.null <- rnorm(500)
  stat.obs <- max(stat.null) + 10

  expect_equal(mean(stat.null >= stat.obs), 0)
  expect_gt(CMRSS:::perm_pvalue(stat.null, stat.obs), 0)
})


test_that("the correction is conservative and shrinks as null.max grows", {
  skip_unless_helper()

  set.seed(20260907)
  stat.null <- rnorm(1000)
  stat.obs <- quantile(stat.null, 0.95, names = FALSE)

  uncorrected <- mean(stat.null >= stat.obs)
  corrected <- CMRSS:::perm_pvalue(stat.null, stat.obs)

  # The correction can only move the p-value up, so a test using it cannot
  # reject where the uncorrected test would not.
  expect_gte(corrected, uncorrected)

  # The gap is bounded by 1 / (m + 1), so the correction matters at small
  # null.max and washes out at large null.max.  This is the reason it is a
  # numerical-output change rather than a cosmetic one: the CMRSS default
  # null.max = 10^4 puts the shift at about 1e-4, but the confidence
  # intervals in the profiling memo ran at null.max = 1000.
  expect_lte(corrected - uncorrected, 1 / (length(stat.null) + 1))
})


test_that("the corrected p-value controls the type I error rate", {
  skip_unless_helper()
  skip_on_cran()

  # The substantive point of the correction, checked rather than asserted.
  # Under the null the observed statistic is exchangeable with the draws,
  # so we simulate that directly: draw m + 1 values from one distribution,
  # call the first the observation and the rest the reference.
  #
  # Phipson and Smyth's claim is that the uncorrected proportion rejects
  # more often than alpha while the corrected one does not.  With m = 19
  # draws and alpha = 0.05 the effect is large enough to see in 4000
  # replicates: the uncorrected rule rejects whenever no draw exceeds the
  # observation, which happens with probability 1/20 = 0.05, and it also
  # rejects on ties, pushing the rate above alpha.

  set.seed(4321)
  m <- 19
  alpha <- 0.05
  reps <- 4000

  reject <- replicate(reps, {
    draws <- rnorm(m + 1)
    obs <- draws[1]
    null <- draws[-1]
    c(uncorrected = mean(null >= obs) <= alpha,
      corrected = CMRSS:::perm_pvalue(null, obs) <= alpha)
  })

  rate_uncorrected <- mean(reject["uncorrected", ])
  rate_corrected <- mean(reject["corrected", ])

  # Monte Carlo standard error of a rate near 0.05 over 4000 replicates is
  # about 0.0034, so three standard errors is about 0.01.
  expect_lte(rate_corrected, alpha + 0.01)
  expect_gt(rate_uncorrected, rate_corrected)
})


test_that("pval_comb_block still uses the uncorrected form", {
  # A tripwire rather than an endorsement.  This pins the current behavior
  # so that wiring in the correction cannot happen silently: whoever makes
  # that change has to come here and delete this test, which is the moment
  # to check that the R&R replication run has already happened.
  # From 0.2.21 the uncorrected share is computed by share_at_least()
  # (R/ties.R), which counts null draws tied with stat.min up to rounding
  # but adds nothing to numerator or denominator.
  src <- deparse(CMRSS:::pval_comb_block)
  expect_true(any(grepl("share_at_least(stat.null, stat.min)", src, fixed = TRUE)))
})
