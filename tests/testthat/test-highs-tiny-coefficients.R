# HiGHS drops constraint-matrix entries with absolute value at most 1e-9 and
# warns once per solve ("LP matrix packed vector contains ... values ...
# ignored"). Stephenson scores with a large s produce such entries: in a site
# of 54 teachers with s = 54, choose(r - 1, 53) / choose(53, 53) is 0 for
# every rank but the top, and after the scores are standardized the near-zero
# ranks of other statistics fall below 1e-9. A confidence interval solves
# hundreds of these problems, so one cmrss() call printed hundreds of
# warnings. Because HiGHS ignores the entries anyway, removing them before
# the call changes no solution; it only removes the warnings.

test_that("HiGHS solves on electric_teachers with large s raise no solver warning", {
  skip_if_not(solver_available("highs"), "HiGHS not available")
  data(electric_teachers, package = "CMRSS", envir = environment())
  d <- electric_teachers
  nb <- as.vector(table(factor(d$Site)))
  ml <- CMRSS:::stephenson_methods(c(2, 11, 59), nb)
  set.seed(1)
  solver_warnings <- character(0)
  withCallingHandlers(
    pval_comb_block(d$TxAny, d$gain, k = 82, c = 0, factor(d$Site), ml,
                    null.max = 200, opt.method = "ILP_highs"),
    warning = function(w) {
      if (!inherits(w, "cmrss_ties_warning")) {
        solver_warnings <<- c(solver_warnings, conditionMessage(w))
      }
      invokeRestart("muffleWarning")
    }
  )
  expect_length(solver_warnings, 0)
})

test_that("dropping entries at or below 1e-9 keeps every other entry", {
  kept <- CMRSS:::drop_tiny_entries(i = 1:4, j = c(1, 1, 2, 2),
                                    x = c(1, 1e-12, -1e-10, -0.5))
  expect_equal(kept$i, c(1L, 4L))
  expect_equal(kept$x, c(1, -0.5))
})
