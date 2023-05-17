#' Generate randomization distribution of the rank score statistic
#'
#'  Generate the null distribution of the given test statistic for an
#'  experiment with m treated out of n units. Same with null_dist in RIQITE.
#'
#' @param Z.perm Matrix of possible treatment assignments. NULL will generate
#' @param score The list of possible scores of the units (under a null,
#'  invariant under treatment assignment). If null will generate from the
#'  passed method.list and given n, m. 
#' @param nperm Number of random permutations
#' @Z.perm User-specified matrix for random permutations

null_dist <- function(n, m, method.list = NULL, score = NULL, 
                      nperm = 10^5, Z.perm = NULL){
  if(is.null(score)){
    score = rank_score( n, method.list )
  }
  
  if(is.null(Z.perm)){
    Z.perm = assign_CRE(n, m, nperm)
  } else {
    stopifnot( is.matrix(Z.perm) )
  }
  
  stat_null = rep(NA, ncol(Z.perm))
  for(iter in 1:ncol(Z.perm)){
    stat_null[iter] = sum( score[ Z.perm[, iter] == 1 ] )
  }
  
  return(stat_null)
}
