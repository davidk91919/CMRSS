# Tests for rank score functions

test_that("rank_score with scale=FALSE produces raw scores", {
  n <- 100

  # Wilcoxon: should be 1:n
  method.wil <- list(name = "Wilcoxon", scale = FALSE)
  expect_equal(rank_score(n, method.wil), c(1:n))

  # Stephenson: should be choose(1:n - 1, s - 1) without normalization.
  # Since 0.2.19 the scores are divided by the top score unless
  # divide_by_top = FALSE (test-stephenson-divide-by-top.R).
  s <- 3
  method.ste <- list(name = "Stephenson", s = s, scale = FALSE, divide_by_top = FALSE)
  expected_ste <- choose(c(1:n) - 1, s - 1)
  expect_equal(rank_score(n, method.ste), expected_ste)
})


test_that("rank_score with scale=TRUE produces standardized scores", {
  n <- 100

  # Wilcoxon with scale=TRUE should be z-score standardized
  method.wil <- list(name = "Wilcoxon", scale = TRUE)
  result <- rank_score(n, method.wil)
  expect_equal(mean(result), 0, tolerance = 1e-10)
  expect_equal(sd(result), 1, tolerance = 1e-10)

  # Stephenson with scale=TRUE should also be standardized
  s <- 4
  method.ste <- list(name = "Stephenson", s = s, scale = TRUE)
  result_ste <- rank_score(n, method.ste)
  expect_equal(mean(result_ste), 0, tolerance = 1e-10)
  expect_equal(sd(result_ste), 1, tolerance = 1e-10)
})


test_that("rank_score_riqite matches RIQITE-style normalization (divide by max)", {
  n <- 50

  # Wilcoxon: rank_score_riqite should divide by max
  method.wil <- list(name = "Wilcoxon")
  result_wil <- rank_score_riqite(n, method.wil)
  expected_wil <- c(1:n) / n
  expect_equal(result_wil, expected_wil)

  # Stephenson: rank_score_riqite should divide by max
  s <- 3
  method.ste <- list(name = "Stephenson", s = s)
  result_ste <- rank_score_riqite(n, method.ste)
  raw_ste <- choose(c(1:n) - 1, s - 1)
  expected_ste <- raw_ste / max(raw_ste)
  expect_equal(result_ste, expected_ste)
})


test_that("Polynomial scores work correctly", {
  n <- 30

  # Polynomial with r=2 and std=FALSE
  method.poly <- list(name = "Polynomial", r = 2, std = FALSE, scale = FALSE)
  result <- rank_score(n, method.poly)
  expected <- (c(1:n))^(2 - 1)  # = 1:n
  expect_equal(result, expected)

  # Polynomial with r=3 and std=TRUE
  method.poly.std <- list(name = "Polynomial", r = 3, std = TRUE, scale = FALSE)
  result_std <- rank_score(n, method.poly.std)
  expected_std <- (c(1:n) / (n + 1))^(3 - 1)
  expect_equal(result_std, expected_std)
})


test_that("rank_score_riqite errors for DIM", {
  expect_error(rank_score_riqite(10, list(name = "DIM")),
               "Can't calculate rank scores for DIM")
})
