#' Inference about quantiles of individual treatment effects from a formula
#'
#' `cmrss()` reads an experiment from a formula and a data frame and returns
#' lower confidence bounds for the sorted individual treatment effects, and,
#' when `quantile` is given, a p-value for the hypothesis that one of those
#' effects is at most `c`. It builds the rank statistics, converts a
#' proportion into the `k` each underlying function expects, and calls
#' [com_conf_quant_larger_cre()] and [comb_p_val_cre()] for a completely
#' randomized experiment or [com_block_conf_quant_larger()] and
#' [pval_comb_block()] for a block-randomized one.
#'
#' @param formula `outcome ~ treatment` for a completely randomized
#'   experiment, or `outcome ~ treatment | block` for a block-randomized one.
#'   The treatment must be coded 0/1 or `TRUE`/`FALSE`. Writing `-outcome`
#'   on the left multiplies every outcome by -1, which turns the
#'   lower bounds into upper bounds on the original effects: the way to ask
#'   whether anyone could have been harmed.
#' @param data A data frame holding the variables in `formula`.
#' @param quantile Optional proportion above 0 and at most 1. When given,
#'   `cmrss()` tests whether the effect at that proportion of the units in
#'   `set` is at most `c`. For example, `quantile = 0.9, set = "treat"` with
#'   164 treated units tests the floor(0.9 x 164) = 147th smallest of their
#'   effects. The `k` used is reported in the result. `quantile = 1` tests
#'   the largest effect; with `-outcome` on the left and `c = 0` it tests the
#'   null hypothesis that no unit in `set` was harmed.
#' @param c The threshold in the tested hypothesis, and the value against
#'   which `print()` counts bounds.
#' @param set Which units' effects to bound: `"treat"`, `"control"`, or
#'   `"all"`.
#' @param scores `"polynomial"` (the default) or `"stephenson"`. A
#'   polynomial score with parameter s gives the unit at rank r, counting
#'   from the smallest outcome, the score (r/(n + 1))^(s - 1). The paper
#'   behind this package writes this parameter zeta. A Stephenson score with
#'   parameter s gives the unit the score choose(r - 1, s - 1). In a
#'   block-randomized experiment both are computed within each block, with
#'   that block's size as n. With the same s the two scores are close
#'   whenever r is much larger than s; see the section "Choosing s".
#' @param s Score parameters, one rank statistic per value; the test
#'   combines them. They are zeta for polynomial scores and s for Stephenson
#'   scores. The default depends on the design and is described in the
#'   section "Choosing s" below. With Stephenson scores, an s above a
#'   block's size is lowered to that size; polynomial scores are never 0 and
#'   need no change.
#' @param alpha One minus the confidence level. Also enters the default `s`.
#' @param nperm Number of simulated random assignments used to approximate
#'   the randomization distribution.
#' @param tol Precision of the confidence bounds.
#' @param opt.method Solver for block-randomized experiments; see
#'   [pval_comb_block()].
#'
#' @section Choosing s:
#'
#' A polynomial rank statistic with parameter s gives the unit at rank r,
#' counting from the smallest outcome, the score (r/(n + 1))^(s - 1), and
#' sums the scores of the treated units. Write u = r/(n + 1). The top 1/s of
#' the ranks carry about the fraction 1 - (1 - 1/s)^s of the total score:
#' 0.75 at s = 2, and no less than 1 - 1/e = 0.63 for any s. So s = 2, whose
#' score is proportional to the rank and gives the Wilcoxon statistic,
#' detects effects shared by most units, and a large s detects large effects
#' confined to a few units. A Stephenson score with the same s,
#' choose(r - 1, s - 1), divided by the top unit's score, is the product of
#' (r - j)/(n - j) over j = 1, ..., s - 1, which is close to (r/n)^(s - 1)
#' when r is much larger than s. So the two families put nearly the same
#' relative weight on each rank in large samples.
#'
#' The default uses three values: 2, s_max, and their geometric middle,
#' sqrt(2 x s_max) rounded to an integer. s_max depends on the number of
#' treated units, m, and the number of units, n, through two steps.
#'
#' First, q_min is the fewest treated units whose effects could be detected
#' at level `alpha`. Suppose only q treated units respond, their effects give
#' them the q largest outcomes, and no other unit has an effect. The data
#' most unlike the null hypothesis of no effect are then those in which the
#' q largest outcomes all belong to treated units. Under no effect that
#' happens with probability
#' (m/n)((m - 1)/(n - 1)) ... ((m - q + 1)/(n - q + 1)). q_min is the
#' smallest q for which this probability is at most `alpha`. With n = 200,
#' m = 100 and alpha = 0.05, q = 4 gives 0.061 and q = 5 gives 0.030, so
#' q_min = 5.
#'
#' Second, s_max = 4 m / q_min, rounded down and never above n / 2. In a
#' simulation of a completely randomized experiment with n = 200 and
#' m = 100, in which a fraction p of the units responded with a large effect
#' and the rest had none, the most powerful single Stephenson statistic had
#' s between about 2/p and 4/p. With q treated responders, p is about q / m,
#' so 4/p is 4 m / q. In the scenario with about 5 treated responders, power
#' peaked at s = 80 = 4 x 100 / 5. Above about s = n / 2 the scores put
#' nearly all their weight on the few top ranks, and power stopped changing.
#'
#' | Design | q_min | s_max | Default s |
#' |---|---|---|---|
#' | n = 200, m = 100 | 5 | 80 | 2, 13, 80 |
#' | `electric_teachers`: n = 233, m = 164 | 9 | 72 | 2, 12, 72 |
#' | n = 30, m = 15 | 4 | 15 | 2, 5, 15 |
#' | n = 100, m = 20 | 2 | 40 | 2, 9, 40 |
#'
#' In a block-randomized experiment q_min and s_max are computed from the
#' total n and the total m, and the same values of s are used in every
#' block. A polynomial score is never 0, so a block smaller than s keeps
#' every unit in the statistic. With a large s almost all of the block's
#' weight sits on its top unit: in a block of 5 with s = 80, the top unit
#' scores (5/6)^79 = 6 x 10^-7 and the next unit (4/6)^79 = 1 x 10^-14. A
#' Stephenson score with s above the block's size would be 0 for every unit
#' in the block, so with `scores = "stephenson"` such an s is lowered to the
#' block's size n_b. The score choose(r - 1, n_b - 1) is then 1 for the top
#' rank and 0 below it, so the block counts toward the statistic exactly
#' when its highest outcome belongs to a treated unit.
#'
#' Why three values. The combined test rejects when the smallest of the
#' separate p-values is at or below a cutoff, the 5th percentile (for
#' `alpha = 0.05`) of that smallest p-value over random assignments under no
#' effect. With one statistic the cutoff is `alpha`. Each added statistic
#' can only lower the cutoff, and that is the only way it can lower power:
#' the test keeps the smallest p-value, so a statistic that carries no
#' information never makes the smallest p-value larger. The cutoff falls
#' most when the added statistic differs from those already present.
#' Statistics with nearby s values reach small p-values on the same
#' assignments; under no effect, the polynomial statistics with s = 6 and
#' s = 10 have correlation 0.96. So the range of s, not the number of
#' values, sets the cutoff.
#'
#' The table gives simulated power at level 0.05 for polynomial scores with
#' n = 200 and m = 100, standard normal control outcomes, and a fraction of
#' units responding with the stated effect, from 2000 replications (Monte
#' Carlo standard error at most 0.011). The last row gives each set's cutoff.
#'
#' | Units responding, effect | s = 2 | 2, 80 | 2, 13, 80 | 2, 4, 8, 16, 32, 64, 80 |
#' |---|---|---|---|---|
#' | 200, 0.35 | 0.79 | 0.73 | 0.72 | 0.72 |
#' | 100, 0.8 | 0.87 | 0.84 | 0.85 | 0.86 |
#' | 40, 1.6 | 0.60 | 0.71 | 0.79 | 0.80 |
#' | 20, 2.5 | 0.32 | 0.73 | 0.77 | 0.78 |
#' | 10, 4 | 0.15 | 0.64 | 0.60 | 0.59 |
#' | 5, 6 | 0.10 | 0.24 | 0.22 | 0.21 |
#' | cutoff | 0.050 | 0.027 | 0.022 | 0.019 |
#'
#' Three values came within 0.02 of the seven-value grid in every scenario.
#' Two values lost 0.79 - 0.71 = 0.08 when 40 units responded, because no
#' value near the best s for that scenario, about 10, was present. Relative
#' to s = 2 alone, the three-value default gave up 0.79 - 0.72 = 0.07 when
#' every unit responded and gained 0.60 - 0.15 = 0.45 when 10 responded.
#' Stephenson scores with the same s values, on the same simulated data,
#' had power within 0.01 of these in every cell.
#'
#' These numbers come from one outcome distribution, effects that are
#' either large or zero, and the null of no effect rather than the quantile
#' hypotheses; they show a pattern, not a guarantee. Supply your own `s`
#' when you have reason to expect effects confined to a smaller or larger
#' share of units than the default covers.
#'
#' @return An object of class `"cmrss"`: a list with `bounds`, a data frame
#'   with one row per sorted effect in `set` (columns `k`, `proportion`,
#'   `lower`), `test` (when `quantile` is given: `k`, `c`, `p.value`), and
#'   the settings used.
#'
#' @examples
#' data(electric_teachers)
#' # Completely randomized analysis, ignoring sites; few permutations so the
#' # example runs quickly. The outcome has many ties, so a warning is given.
#' fit <- cmrss(gain ~ TxAny, data = electric_teachers, set = "treat",
#'              nperm = 500, tol = 0.1)
#' fit
#' head(fit$bounds)
#'
#' # Could anyone have been harmed? Effects on -gain are minus the effects
#' # on gain, so a lower bound above 0 here means some teachers lost.
#' harm <- cmrss(-gain ~ TxAny, data = electric_teachers, set = "all",
#'               nperm = 500, tol = 0.1)
#' sum(harm$bounds$lower > 0)
#'
#' \donttest{
#' # Block-randomized by site (needs the highs package)
#' if (requireNamespace("highs", quietly = TRUE)) {
#'   cmrss(gain ~ TxAny | Site, data = electric_teachers, set = "treat",
#'         nperm = 500, opt.method = "ILP_highs")
#' }
#' }
#'
#' @seealso [com_conf_quant_larger_cre()], [com_block_conf_quant_larger()]
#' @export
cmrss <- function(formula, data, quantile = NULL, c = 0, set = "all",
                  scores = c("polynomial", "stephenson"), s = NULL,
                  alpha = 0.05, nperm = 10^4, tol = 0.01,
                  opt.method = "ILP_auto") {
  set <- match.arg(set, c("all", "treat", "control"))
  scores <- match.arg(scores)
  d <- parse_cmrss_formula(formula, data)
  y <- d$y; z <- d$z; block <- d$block
  n <- length(y); m <- sum(z)

  if (!is.null(quantile) &&
      (!is.numeric(quantile) || length(quantile) != 1L ||
       quantile <= 0 || quantile > 1)) {
    stop("quantile must be a single proportion above 0 and at most 1.")
  }

  # One warning per call: check here, then muffle the same warning from the
  # functions cmrss() calls.
  check_outcome_ties(y)
  quietly <- function(expr) {
    withCallingHandlers(expr, cmrss_ties_warning = function(w) {
      invokeRestart("muffleWarning")
    })
  }

  if (is.null(s)) s <- default_score_parameters(n, m, alpha)
  nb <- if (is.null(block)) n else as.vector(table(block))
  ml_all <- if (scores == "polynomial") polynomial_methods(s, nb) else
    stephenson_methods(s, nb)
  # The completely randomized functions take one specification per statistic.
  ml <- if (is.null(block)) lapply(ml_all, `[[`, 1) else ml_all

  # The test comes before the bounds so that a caller who sets a seed gets
  # the same p-value as a direct call to the underlying function.
  test <- NULL
  if (!is.null(quantile)) {
    test <- cmrss_test(y, z, block, quantile, c, set, ml, nperm, opt.method,
                       quietly)
  }

  lower <- quietly(if (is.null(block)) {
    com_conf_quant_larger_cre(z, y, ml, nperm = nperm, set = set,
                              alpha = alpha, tol = tol)
  } else {
    com_block_conf_quant_larger(z, y, block, set = set,
                                methods.list.all = ml,
                                opt.method = opt.method, null.max = nperm,
                                tol = tol, alpha = alpha)
  })
  size <- length(lower)

  structure(list(
    bounds = data.frame(k = seq_len(size), proportion = seq_len(size) / size,
                        lower = lower),
    test = test, set = set, scores = scores, s = s, c = c, alpha = alpha,
    nperm = nperm, n = n, m = m, blocked = !is.null(block),
    outcome = paste(deparse(formula[[2]]), collapse = " "), call = match.call()
  ), class = "cmrss")
}

# Read outcome, treatment and optional block from `y ~ z` or `y ~ z | b`.
parse_cmrss_formula <- function(formula, data) {
  env <- environment(formula)
  rhs <- formula[[3]]
  block_expr <- NULL
  if (is.call(rhs) && identical(rhs[[1]], as.name("|"))) {
    block_expr <- rhs[[3]]
    rhs <- rhs[[2]]
  }
  y <- eval(formula[[2]], data, env)
  z <- eval(rhs, data, env)
  block <- if (is.null(block_expr)) NULL else eval(block_expr, data, env)

  missing <- is.na(y) | is.na(z)
  if (!is.null(block)) missing <- missing | is.na(block)
  if (any(missing)) {
    stop(sprintf(paste0(
      "%d rows have missing values in the outcome, treatment or block. ",
      "Remove them, or impute them, before calling cmrss()."), sum(missing)))
  }
  if (is.logical(z)) z <- as.numeric(z)
  if (!is.numeric(z) || !all(z %in% c(0, 1))) {
    stop("The treatment must be coded 0 and 1, or TRUE and FALSE.")
  }
  if (!is.numeric(y)) stop("The outcome must be numeric.")
  list(y = as.numeric(y), z = as.numeric(z),
       block = if (is.null(block)) NULL else factor(block))
}

# p-value for H: the k-th smallest effect in `set` is at most c, with k set
# by `quantile`. Each underlying function counts k its own way:
# comb_p_val_cre over all n units, pval_comb_block over the treated only.
cmrss_test <- function(y, z, block, quantile, c, set, ml, nperm, opt.method,
                       quietly) {
  n <- length(y)
  # One side: effects among the units coded 1 in zz, with k counted among
  # them (k_treated) or, when that k is not a treated quantile, p = 1.
  one_side <- function(zz, yy, k_treated) {
    m <- sum(zz)
    if (k_treated < 1) return(1)
    quietly(if (is.null(block)) {
      comb_p_val_cre(zz, yy, k = n - m + k_treated, c = c, ml, nperm = nperm)
    } else {
      pval_comb_block(zz, yy, k = k_treated, c = c, block, ml,
                      null.max = nperm, opt.method = opt.method,
                      statistic = FALSE)
    })
  }
  size <- switch(set, treat = sum(z), control = n - sum(z), all = n)
  k <- floor(quantile * size)
  if (k < 1) {
    stop(sprintf("quantile = %s of %d units is below the smallest effect.",
                 format(quantile), size))
  }
  p <- switch(set,
    treat = one_side(z, y, k),
    control = one_side(1 - z, -y, k),
    # The k-th smallest of all n effects is at most c only if at most n - k
    # treated and at most n - k control effects exceed c. Each side is tested
    # and the smaller p-value doubled, the Bonferroni step that
    # com_conf_quant_larger_cre and com_block_conf_quant_larger use for
    # set = "all".
    all = min(1, 2 * min(one_side(z, y, k - (n - sum(z))),
                         one_side(1 - z, -y, k - sum(z))))
  )
  list(k = k, c = c, p.value = p)
}
