#' rank_score
#'
#' @description Calculating the scores of n units (ordered from low to high score) given specified test statistic, not allowing for ties.
#'
#' @param n Sample size
#' @param method.list List of method. Standarizing and scaling is available for Polynomial rank score. Scaling is available for Wilcoxon and Stephenson.

rank_score   <- function(n,
                         method.list = list(
                           name = "Polynomial", r, std = TRUE, scale = FALSE)
){
  if(method.list$name == "Polynomial"){
    r = method.list$r
    if(method.list$std == TRUE) {
      score = (c(1:n) / n)^(r-1)
    } else {
      score = (c(1:n))^(r-1)
    }
    if(method.list$scale == TRUE){
      score = scale(score)
    }
  }
  if(method.list$name == "Stephenson"){
    score = choose( c(1:n) - 1, method.list$s - 1)
    if(method.list$scale == TRUE){
      score = scale(score)
    }
  }
  if(method.list$name == "Wilcoxon"){
    score = c(1:n)
    if(method.list$scale == TRUE){
      score = scale(score)
    }
  }
  score = score / max(score)
  return(score)
}
