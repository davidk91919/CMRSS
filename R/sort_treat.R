#' sort_treat
#'
#' @noRd





sort_treat <- function(Y, Z){
  r = rank(Y, ties.method = "first")
  ind.sort = sort.int(r, index.return = TRUE)$ix
  ind.sort.treat = ind.sort[Z[ind.sort] == 1]
  return(ind.sort.treat)
}
