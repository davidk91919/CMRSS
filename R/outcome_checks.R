# Warning about binary and heavily tied outcomes.
#
# Every rank in the package uses rank(..., ties.method = "first"), so a tie
# between a treated and a control unit is broken by which row comes first in
# the data, and reordering rows changes results (issue #5). Changing the
# rule would move the numbers in the combined_stephenson_tests paper, so
# until that paper is published the user-facing functions warn instead.
# Rule chosen by Jake on 2026-10-04: two distinct values, or one value shared
# by more than 5 percent of the units. A value held by one unit is not a tie,
# however small the sample.

#' Warn when ties in the outcome make results depend on row order
#'
#' @param Y Outcome vector.
#' @param threshold Largest share of units that may share one value before
#'   the warning is given.
#'
#' @return `Y`, invisibly. Called for its warning, which has class
#'   `"cmrss_ties_warning"`.
#' @keywords internal
check_outcome_ties <- function(Y, threshold = 0.05) {
  Y <- Y[!is.na(Y)]
  n <- length(Y)
  if (n == 0) return(invisible(Y))
  counts <- table(Y)
  why <- NULL
  if (length(counts) == 2) {
    why <- "The outcome takes only two values."
  } else if (max(counts) >= 2 && max(counts) / n > threshold) {
    why <- sprintf("%d of %d units (%.0f percent) share the outcome value %s.",
                   max(counts), n, 100 * max(counts) / n,
                   names(counts)[which.max(counts)])
  }
  if (!is.null(why)) {
    msg <- paste(
      why,
      "Tied outcomes are ranked by their order in the data, so reordering the",
      "rows can change p-values and confidence bounds (see",
      "https://github.com/davidk91919/CMRSS/issues/5).",
      if (length(counts) == 2)
        "Rank-based quantile inference is not designed for binary outcomes."
    )
    warning(structure(class = c("cmrss_ties_warning", "warning", "condition"),
                      list(message = msg, call = sys.call(-1))))
  }
  invisible(Y)
}
