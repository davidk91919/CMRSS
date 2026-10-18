# In the program HiGHS_sol_com() writes for a stratified p-value, one
# variable holds each statistic being combined (eta_h) and one the largest
# of the standardized statistics (theta), which the program minimizes.
# Through CMRSS 0.2.20 those H + 1 variables had no bounds, and HiGHS
# 1.14.0.2 sometimes stopped on such a program with a C++ error,
# "vector::reserve" on Linux and "vector" on macOS. In a simulation in 20
# strata of 50 units on the Illinois Campus Cluster
# (combined_stephenson_tests, code/simulations/zeta_regret/README.md,
# section 8.8) one program failed with presolve on and again with presolve
# off, so the fallback of 0.2.20 could not rescue it. The same program
# solves with either setting once those variables have finite bounds.
#
# The bounds come from the program itself. The program sets one x to 1 in
# each stratum (in the relaxation, the x's of a stratum are weights that sum
# to 1), so eta_h is a sum over strata of one coefficient per stratum, or of
# a weighted average of them. eta_h therefore lies between the sum of the
# strata's smallest coefficients and the sum of their largest, and theta,
# the largest (eta_h - mu_h) / sigma_h, lies between the largest of those
# standardized lower ends and the largest of the upper ends. No choice of
# the x's is excluded, so the optimum cannot change. The tests check this
# against every choice in a small program, against the optimum of the same
# programs without bounds, and on the program that failed.

## A stratified experiment small enough to list every choice the program
## can make: 4 strata of 8 units with 4 treated in each, so 5 options per
## stratum and 5^4 = 625 choices, and two statistics combined.
bounds_sre <- function() {
  set.seed(20261018)
  block <- factor(rep(1:4, each = 8))
  Z <- unlist(lapply(1:4, function(b) sample(rep(0:1, each = 4))))
  Y <- rnorm(32) + Z
  nb <- as.vector(table(block))
  list(Z = Z, Y = Y, block = block, ml = CMRSS:::polynomial_methods(c(2, 6), nb))
}

## The arguments pval_comb_block() (comb.method = 1) passes to
## HiGHS_sol_com() for its k, built the way pval_comb_block() builds them.
sol_com_args <- function(d, k, exact = TRUE, c = 0) {
  bs <- CMRSS:::summary_block(d$Z, d$block)
  weight <- CMRSS:::weight_scheme(bs, "asymp.opt")
  scores <- lapply(d$ml, function(ml) CMRSS:::score_all_blocks(bs$nb, ml))
  list(Z = d$Z, block = bs$block, weight = weight,
       coeflists = CMRSS:::comb_matrix_block(d$Z, d$Y, bs$block, c, d$ml, scores, bs),
       p = sum(d$Z) - k,
       ms_list = CMRSS:::mu_sigma_list(d$Z, bs$block, weight, d$ml, scores, bs),
       exact = exact, block.sum = bs)
}

test_that("the bounds contain every value the statistics can take", {
  d <- bounds_sre()
  a <- sol_com_args(d, k = 8)
  H <- length(a$coeflists); B <- length(a$coeflists[[1]]); mb <- a$block.sum$mb
  ## the coefficient of each option of stratum b in statistic h's row of
  ## the program (HiGHS_sol_com(), constraint 3)
  coef <- function(h, b) a$weight[b] / mb[b] * a$coeflists[[h]][[b]][2, ]
  ## every choice of one option per stratum, and the statistics it gives
  choices <- as.matrix(expand.grid(lapply(seq_len(B), function(b) seq_along(coef(1, b)))))
  eta <- sapply(seq_len(H), function(h) {
    apply(choices, 1, function(j) sum(vapply(seq_len(B), function(b) coef(h, b)[j[b]], 0)))
  })
  theta <- apply(sweep(sweep(eta, 2, a$ms_list$mu), 2, a$ms_list$sigma, "/"), 1, max)
  values <- cbind(eta, theta)

  bnd <- CMRSS:::highs_com_bounds(a$coeflists, a$weight, mb, a$ms_list)
  expect_true(all(is.finite(c(bnd$lower, bnd$upper))))
  ## no value lies outside its bounds
  expect_true(all(sweep(values, 2, bnd$lower, ">=")))
  expect_true(all(sweep(values, 2, bnd$upper, "<=")))
  ## and each eta's bounds are its smallest and largest values, as is
  ## theta's upper bound, up to the padding: a millionth of the value
  ## against rounding, plus 0.001 so that a variable whose smallest and
  ## largest values coincide is not nearly fixed (see the test of equal
  ## coefficients below). theta's lower bound can sit below its smallest
  ## value, because the statistics can reach their smallest values at
  ## different choices.
  pad <- function(x) 1e-3 + 1e-6 * pmax(1, abs(x))
  lo_eta <- apply(eta, 2, min)
  hi_all <- unname(apply(values, 2, max))
  expect_equal(bnd$lower[1:H], lo_eta - pad(lo_eta), tolerance = 1e-12)
  expect_equal(bnd$upper, hi_all + pad(hi_all), tolerance = 1e-12)
})

test_that("bounding the statistics changes no optimum, exact or relaxed", {
  skip_if_not(solver_available("highs"), "HiGHS not available")
  d <- bounds_sre()
  solve_at <- function(k, exact) do.call(CMRSS:::HiGHS_sol_com, sol_com_args(d, k, exact))$obj
  grid <- expand.grid(k = c(1, 4, 8, 12, 16), exact = c(TRUE, FALSE))
  bounded <- mapply(solve_at, grid$k, grid$exact)
  ## the same programs with the bounds removed, as CMRSS 0.2.20 wrote them
  local_mocked_bindings(highs_com_bounds = function(coeflists, ...) {
    H <- length(coeflists)
    list(lower = rep(-Inf, H + 1), upper = rep(Inf, H + 1))
  })
  unbounded <- mapply(solve_at, grid$k, grid$exact)
  expect_equal(bounded, unbounded, tolerance = 1e-9)
})

test_that("a program whose strata's coefficients are all equal solves to its optimum, presolve on and off", {
  skip_if_not(solver_available("highs"), "HiGHS not available")
  ## With c = 100 every treated unit's Y - c ranks below every control, so
  ## the treated units take the lowest ranks in their stratum whatever the
  ## number exempted, and each stratum's coefficients are one number
  ## repeated. Each eta's smallest and largest values then coincide. With
  ## bounds widened only by a millionth of that value, HiGHS 1.14.0.2
  ## returned the status "Solve error" and objective 0 on such a program
  ## with presolve on, and "Infeasible" with it off.
  d <- bounds_sre()
  a <- sol_com_args(d, k = sum(d$Z), c = 100)
  mb <- a$block.sum$mb
  ## (equal up to rounding: each is a sum of the same scores in a different order)
  expect_true(all(vapply(a$coeflists, function(cl)
    all(vapply(cl, function(m) diff(range(m[2, ])) < 1e-12, TRUE)), TRUE)))
  ## the only choice: every stratum exempts no unit (p = 0)
  eta <- vapply(a$coeflists, function(cl)
    sum(vapply(seq_along(cl), function(b) a$weight[b] / mb[b] * cl[[b]][2, 1], 0)), 0)
  optimum <- max((eta - a$ms_list$mu) / a$ms_list$sigma)

  real <- highs::highs_solve
  for (presolve in c("on", "off")) {
    local_mocked_bindings(highs_solve = function(..., control) {
      control$presolve <- presolve
      real(..., control = control)
    }, .package = "highs")
    expect_no_warning(result <- do.call(CMRSS:::HiGHS_sol_com, a))
    expect_equal(result$obj, optimum, tolerance = 1e-9, label = paste("presolve", presolve))
  }
})

test_that("the program that stopped HiGHS with presolve on and off solves once bounded", {
  skip_if_not(solver_available("highs"), "HiGHS not available")
  ## The arguments of HiGHS_sol_com() in replication 5541 of the cell in 20
  ## strata of 50 where the top fifth on a covariate gains: k = 450, five
  ## polynomial degrees, saved on the cluster. There HiGHS stopped on the
  ## unbounded program under every setting tried; on macOS it stops with
  ## presolve off. HiGHS is made to run with presolve off here, so that
  ## without the bounds the test fails on macOS too.
  a <- readRDS(test_path("fixtures", "highs_unbounded_failure.rds"))
  real <- highs::highs_solve
  local_mocked_bindings(highs_solve = function(..., control) {
    control$presolve <- "off"
    real(..., control = control)
  }, .package = "highs")
  expect_no_warning(result <- do.call(CMRSS:::HiGHS_sol_com, a))
  ## the optimum HiGHS reaches on the unbounded program with presolve on
  expect_equal(result$obj, 1.19090017031, tolerance = 1e-9)
})
