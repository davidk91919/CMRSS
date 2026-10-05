#' Print a cmrss result
#'
#' Reports how many of the sorted effects have a lower confidence bound above
#' `c`: if the bound for the k-th smallest of N effects is above c, so are
#' the bounds for the k-th through N-th, and at least N - k + 1 units have
#' effects above c.
#'
#' @param x An object returned by [cmrss()].
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.cmrss <- function(x, ...) {
  group <- switch(x$set, treat = "treated", control = "control", all = "")
  size <- nrow(x$bounds)
  n_above <- sum(x$bounds$lower > x$c)
  cat(sprintf("cmrss: %s experiment, %d units, %d treated\n",
              if (x$blocked) "block-randomized" else "completely randomized",
              x$n, x$m))
  cat(sprintf("%s scores with parameters %s; %d simulated assignments\n",
              if (x$scores == "stephenson") "Stephenson" else "Polynomial",
              paste(x$s, collapse = ", "), x$nperm))
  # Naming the outcome matters when it is -y: the count is then of units
  # whose effect on y is below -c.
  cat(sprintf(paste0("With %s percent confidence, at least %d of %d %s",
                     "units have effects on %s above %s.\n"),
              format(100 * (1 - x$alpha)), n_above, size,
              if (nzchar(group)) paste0(group, " ") else "", x$outcome,
              format(x$c)))
  if (!is.null(x$test)) {
    cat(sprintf(paste0("Test that effect number k = %d, counting up from ",
                       "the smallest of the %d %seffects, is at most %s: ",
                       "p = %s\n"),
                x$test$k, size, if (nzchar(group)) paste0(group, " ") else "",
                format(x$test$c), format(signif(x$test$p.value, 3))))
  }
  invisible(x)
}
