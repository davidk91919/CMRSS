################################################################################
# Letting a caller cap Gurobi's thread count
#
# Why this exists. Gurobi's default is to use every core on the machine for a
# single solve. That is right for one solve at a time and wrong when the caller
# is already parallel: a foreach loop over 14 workers, each solving with 14
# threads, puts about 150 runnable threads on 14 cores, and the workers spend
# their time queueing behind each other rather than solving. Measured on an
# Apple M3 Max running the combined_stephenson_tests SRE simulation, the load
# average reached 158 on a 14-core machine and the workers together drew 581
# percent of CPU where 14 busy cores would be 1400 percent.
#
# The caller is the only one who knows whether it is already parallel, so the
# thread count belongs to the caller. The option CMRSS.gurobi.threads carries
# it. Unset, which is the default, nothing changes and Gurobi decides as it
# always did.
#
# What must stay true is that this is a resource setting and not a modelling
# one: capping threads must not change the solution. The last test here is the
# one that matters, and it is skipped when Gurobi is not installed.
################################################################################

## Set the option for the duration of one expression and restore it after,
## without taking a dependency on withr, which CMRSS does not declare.
with_threads <- function(value, code) {
  old <- options(CMRSS.gurobi.threads = value)
  on.exit(options(old), add = TRUE)
  force(code)
}

## The internal solvers are where the params list is built.
gurobi_params_source <- function(fn) {
  paste(deparse(body(utils::getFromNamespace(fn, "CMRSS"))), collapse = "\n")
}

solver_fns <- c("Gurobi_sol_com", "Gurobi_sol_stratum_com")


test_that("both Gurobi solvers build their params through the helper", {
  ## Neither solver should assemble a params list of its own, or the option
  ## would reach one call site and not the other.
  for (fn in solver_fns) {
    src <- gurobi_params_source(fn)
    expect_true(grepl("add_thread_param", src, fixed = TRUE),
                info = paste(fn, "does not call add_thread_param, so a parallel",
                             "caller cannot stop it from taking every core"))
  }
})


test_that("the helper is the one place that reads the option", {
  src <- paste(deparse(body(utils::getFromNamespace("add_thread_param", "CMRSS"))),
               collapse = "\n")
  expect_true(grepl("CMRSS.gurobi.threads", src, fixed = TRUE))
})


test_that("the option is unset by default, so behaviour does not change", {
  expect_null(getOption("CMRSS.gurobi.threads"))
})


test_that("the helper turns the option into a params entry", {
  add_threads <- utils::getFromNamespace("add_thread_param", "CMRSS")

  ## Unset: params come back untouched, so Gurobi decides.
  with_threads(NULL, {
    p <- add_threads(list(OutputFlag = 0))
    expect_identical(p, list(OutputFlag = 0))
  })

  ## Set: the count is passed through as an integer.
  with_threads(1, {
    p <- add_threads(list(OutputFlag = 0))
    expect_identical(p$Threads, 1L)
    expect_identical(p$OutputFlag, 0)
  })

  with_threads(4, expect_identical(add_threads(list())$Threads, 4L))

  ## Nonsense is refused rather than passed to the solver, where it would
  ## either be ignored or error somewhere less informative.
  with_threads(0, expect_error(add_threads(list()), "positive"))
  with_threads("many", expect_error(add_threads(list()), "positive"))
})


test_that("capping threads does not change the solution", {
  ## The statistical point. Threads is a resource setting: the same model must
  ## give the same optimum however many cores work on it. If this ever fails,
  ## the cap is not safe to use and the simulation results computed under it
  ## would not be comparable with results computed without it.
  skip_if_not_installed("gurobi")

  set.seed(4)
  n <- 40L
  B <- 4L
  block <- factor(rep(seq_len(B), each = n / B))
  Z <- unlist(lapply(split(seq_len(n), block), function(idx) {
    z <- rep(0, length(idx)); z[sample(length(idx), length(idx) %/% 2)] <- 1; z
  }))
  Y <- round(stats::rnorm(n), 3)

  ## The function draws its own permutations, so without a fixed seed two calls
  ## differ by Monte Carlo noise and the comparison would say nothing about
  ## threads. Seeding first makes the two runs differ only in the thread count.
  run <- function(threads) {
    set.seed(2024)
    with_threads(threads,
      CMRSS::com_block_conf_quant_larger(
        Z = Z, Y = Y, block = block, set = "treat",
        methods.list.all = list(list(
          list(name = "Polynomial", r = 2, std = TRUE, scale = FALSE),
          list(name = "Polynomial", r = 2, std = TRUE, scale = FALSE),
          list(name = "Polynomial", r = 2, std = TRUE, scale = FALSE),
          list(name = "Polynomial", r = 2, std = TRUE, scale = FALSE))),
        weight.name = "asymp.opt", comb.method = 2,
        null.max = 1000, tol = 0.01, alpha = 0.1))
  }

  expect_equal(run(1), run(NULL))
})
