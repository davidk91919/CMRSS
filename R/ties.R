# Comparisons that count values equal up to rounding as tied.
#
# A p-value from CMRSS is the share of null draws at or above the smallest
# value of the statistic under the hypothesis. The two sides of that
# comparison are computed by different routes: the smallest value by a sum
# over the observed ranks or by an integer program, a null draw by a sum
# over a permuted assignment. Two values that are mathematically equal can
# then differ in their last digits, and through 0.2.20 a plain `>=` decided
# the tie by rounding. When the null draw came out lower it was not counted
# and the p-value was too small. The same happens wherever one share is
# compared with another computed differently (1 - 997/1000 is not 3/1000).
#
# Every such comparison in the package goes through the functions here, so
# that the observed value and the null draws are compared under one rule.
# Jake chose the margin on 2026-10-18: two numbers count as equal when they
# differ by less than 1e-9 * max(1, abs(observed number)). Rounding moves
# these values by about 1e-15 of their size, and in the experiments the
# package is used for, values that truly differ are much more than 1e-9
# apart (for shares, at least 1 / nperm).

#' Margin within which two numbers count as tied
#'
#' @param x The observed number (or numbers) the comparison is made against.
#'
#' @return \code{1e-9 * pmax(1, abs(x))}.
#'
#' @keywords internal
tie_tolerance <- function(x) {
  1e-9 * pmax(1, abs(x))
}

#' Share of draws at or above a value, counting near-ties as ties
#'
#' @param draws Numeric vector of null draws.
#' @param value A single observed value.
#'
#' @return The share of \code{draws} that are at least
#'   \code{value - tie_tolerance(value)}.
#'
#' @keywords internal
share_at_least <- function(draws, value) {
  # sum() / length() rather than mean(): mean() of a logical accumulates in
  # long double on x86, so two equal counts need not give bit-identical
  # shares, and the shares are compared with each other later.
  sum(draws >= value - tie_tolerance(value)) / length(draws)
}

#' Share of draws at or below a value, counting near-ties as ties
#'
#' @inheritParams share_at_least
#'
#' @return The share of \code{draws} that are at most
#'   \code{value + tie_tolerance(value)}.
#'
#' @keywords internal
share_at_most <- function(draws, value) {
  sum(draws <= value + tie_tolerance(value)) / length(draws)
}
