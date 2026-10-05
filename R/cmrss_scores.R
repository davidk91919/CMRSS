# Rank-score specifications for cmrss(): the default parameter values and
# the nested lists the stratified functions expect. The reasoning behind the
# default is in the section "Choosing s" of ?cmrss.

#' Fewest treated units whose effects could be detected at level alpha
#'
#' The smallest q for which the probability that the q largest outcomes all
#' belong to treated units, under no effect, is at most `alpha`:
#' (m/n)((m - 1)/(n - 1)) ... ((m - q + 1)/(n - q + 1)) <= alpha.
#'
#' @param n Number of units.
#' @param m Number of treated units.
#' @param alpha Level of the test.
#' @return An integer, or `NA` when no q up to m reaches `alpha`.
#' @keywords internal
min_detectable_treated <- function(n, m, alpha = 0.05) {
  probs <- cumprod((m - seq_len(m) + 1) / (n - seq_len(m) + 1))
  q <- which(probs <= alpha)
  if (length(q) == 0) NA_integer_ else q[1]
}

#' Default score parameters for cmrss(), polynomial zeta or Stephenson s
#'
#' Three values: 2, s_max = min(floor(4 m / q_min), floor(n / 2)), and their
#' geometric middle round(sqrt(2 s_max)).
#'
#' @inheritParams min_detectable_treated
#' @return A sorted integer vector of distinct values.
#' @keywords internal
default_score_parameters <- function(n, m, alpha = 0.05) {
  q_min <- min_detectable_treated(n, m, alpha)
  if (is.na(q_min)) {
    stop(sprintf(paste0(
      "With %d treated units of %d, no rank test can reach alpha = %s: even ",
      "the treated units holding all the top ranks is not that unlikely ",
      "under no effect. Supply s yourself only if you accept that."),
      m, n, format(alpha)))
  }
  s_max <- max(2, min(floor(4 * m / q_min), floor(n / 2)))
  sort(unique(c(2, round(sqrt(2 * s_max)), s_max)))
}

#' Stephenson score lists, one per statistic and block
#'
#' An s larger than a block's size is lowered to that size. The score
#' choose(r - 1, n_b - 1) is then 1 for the block's top rank and 0 below it,
#' so the block still counts instead of scoring 0 everywhere.
#'
#' @param s Stephenson parameters.
#' @param nb Block sizes (one value for a completely randomized experiment).
#' @return A list with one element per value of `s`, each a list with one
#'   score specification per block, the form `methods.list.all` takes.
#' @keywords internal
stephenson_methods <- function(s, nb) {
  lapply(s, function(sv) {
    lapply(nb, function(n_b) list(name = "Stephenson", s = min(sv, n_b),
                                  scale = FALSE))
  })
}

#' Polynomial score lists, one per statistic and block
#'
#' A polynomial score (r / (n_b + 1))^(zeta - 1) is never 0, so zeta is used
#' unchanged in every block.
#'
#' @param zeta Polynomial parameters.
#' @inheritParams stephenson_methods
#' @return A list in the form `methods.list.all` takes.
#' @keywords internal
polynomial_methods <- function(zeta, nb) {
  lapply(zeta, function(z) {
    lapply(nb, function(n_b) list(name = "Polynomial", r = z, std = TRUE,
                                  scale = FALSE))
  })
}
