#' min_stat
#'
#' 
#'
#' @param Z Given assignment vector
#' @param Y Outcome vector
#' @param k k specified as kth largest treatment effect in the null hypothesis
#' @param c c specified as upper bound of kth largest treatment effect in the null hypothesis
#' @param method.list list of method.
#' @param score user-specified score vector, if available
#' @param ind.sort.treat user-specified alignment in treatment vector, if available
#' @noRd

min_stat <- function(Z, Y, k, c, method.list = NULL,
                     score = NULL,
                     ind.sort.treat = NULL){
  n = length(Y)
  m = sum(Z)

  if(is.null(method.list)){   #added
    score= scale(score)
  }

  if(is.null(score)){
    score = rank_score(n, method.list)
  }

  # sort the treated units
  if(is.null(ind.sort.treat)){
    ind.sort.treat = sort_treat(Y, Z)
  }

  # get xi vector
  xi = rep(c, n)
  if(k < n){
    xi[ ind.sort.treat[ ( m + 1 - min(m, n-k) ):m ] ] = Inf
  }

  stat.min = sum( score[ rank( Y - Z * xi, ties.method = "first" )[Z==1] ] )
  stat.min

  return(stat.min)
}
