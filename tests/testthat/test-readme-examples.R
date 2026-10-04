# Tests that the README's worked examples run, and that the error a user
# meets for an out-of-range k in pval_comb_block tells them what k counts.
#
# Why this file exists:
#
# Since 0.2.7, pval_comb_block counts k over the treated units only, so k
# must lie in 1..sum(Z). The README's Quick Start and Example 1 still passed
# k = floor(0.9 * N), counted over all N units, and both stopped with an
# error. Nothing ran the README, so nobody noticed. The error itself cited a
# source line ("R/CMRSS_SRE.R:1034") that a user of the installed package
# cannot see and did not say that k counts treated units.
#
# The README is not installed with the package, so the example test runs
# only from a source checkout (devtools::test()) and skips under R CMD check.

readme_chunks <- function(path) {
  lines <- readLines(path, warn = FALSE)
  starts <- which(lines == "```r")
  ends <- which(lines == "```")
  lapply(starts, function(s) {
    e <- min(ends[ends > s])
    lines[(s + 1):(e - 1)]
  })
}

test_that("every R chunk in README.md runs without error", {
  skip_if_not(solver_available("highs"), "HiGHS not available")
  readme <- test_path("..", "..", "README.md")
  skip_if_not(file.exists(readme), "README.md not available (installed package)")

  chunks <- readme_chunks(readme)
  expect_gt(length(chunks), 0)

  # One environment for the whole README, as a reader running it top to
  # bottom would have. Installation chunks and help lookups are not examples.
  env <- new.env(parent = globalenv())
  for (chunk in chunks) {
    code <- paste(chunk, collapse = "\n")
    if (grepl("install", code) || grepl("^\\s*\\?", code)) next
    if (grepl("ILP_gurobi", code) && !solver_available("gurobi")) next
    err <- tryCatch({
      utils::capture.output(suppressMessages(eval(parse(text = code), envir = env)))
      NULL
    }, error = function(e) conditionMessage(e))
    expect(is.null(err),
           paste0("README chunk failed with: ", err, "\n--- chunk ---\n", code))
  }
})

make_small_sre <- function() {
  set.seed(7)
  s <- 3; n_per <- 8; m_per <- 4
  block <- factor(rep(1:s, each = n_per))
  Z <- rep(0, s * n_per)
  for (i in 1:s) Z[sample(which(block == i), m_per)] <- 1
  Y <- rnorm(s * n_per) + Z
  ml <- list(lapply(1:s, function(i) list(name = "Wilcoxon", scale = FALSE)))
  list(Z = Z, Y = Y, block = block, ml = ml)
}

test_that("out-of-range k: the error says k counts treated units and gives the range", {
  skip_if_not(solver_available("highs"), "HiGHS not available")
  d <- make_small_sre()  # 24 units, 12 treated

  err <- tryCatch(
    pval_comb_block(d$Z, d$Y, k = 21, c = 0, d$block, d$ml,
                    opt.method = "ILP_highs", null.max = 100),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "treated")
  expect_match(err, "between 1 and sum\\(Z\\) = 12")
  # A source line number means nothing to a user of the installed package.
  expect_no_match(err, "R/CMRSS_SRE\\.R")
})

test_that("a k that is not a whole number is refused", {
  skip_if_not(solver_available("highs"), "HiGHS not available")
  d <- make_small_sre()
  # tau_(2.5) is not a quantile; returning a p-value would answer a
  # question nobody can state.
  expect_error(
    pval_comb_block(d$Z, d$Y, k = 2.5, c = 0, d$block, d$ml,
                    opt.method = "ILP_highs", null.max = 100),
    "whole number"
  )
})
