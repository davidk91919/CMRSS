#' @title Combined minimum stat values in CRE
#' 
#' @description Function calculates combined minimum stat values for given quantile null hypothesis H_k,c.
#'
#' @param methods.list lists of methods. Should be list of lists
#' @param combine Function will be used for combining values. Default option is maximum function.
#'
#' @export


combine_min_stat <- function(Z, Y, k, c, methods.list = NULL,
                             ind.sort.treat = NULL,
                             combine = max) {

  if(is.null(scores)){
    K = length(methods.list)
  }else{
    K = length(scores)
  }

  temp <- rep(0, K)

  if(is.null(ind.sort.treat)){
    ind.sort.treat = sort_treat(Y, Z)
  }

  for (i in 1 : K) {
    temp[i] <- min_stat(Z, Y, k, c,
                        method.list = methods.list[[i]],
                        score = scores[[i]],
                        ind.sort.treat = ind.sort.treat)
  }
  result <- combine(temp)

  return(result)
}
