#' @title Calculating all possible test statistic values for every strata
#' @description Generate a list of lists, where elements are all possible test statistic value for every possible null hypothesis in each stratum,
#' using user-specified score function.
#'
#' @param Z n-dimensional treatment assignment vector, for every units.
#' @param Y n-dimensional observed outcome vector, for every units.
#' @param block n-dimensional vector specifying stratum status of each units.
#' @param c A scalar specifying the bounded null hypothesis.
#' @param method.list.all list of every stratified rank sum statistics, where each element(list)s is method for every stratum.List of stratified rank sum statistic. The structure should be list of lists of lists.


###################################################
# function for calculate all possible t_s^h(z, y - z \circ xi)
# by choice of l_s
# note, now score function is not neither Stephenson nor Wilcoxon.
###################################################

comb_matrix_block = function(Z, Y, block, c, methods.list.all){
  if(!is.factor(block)){
    block = as.factor(block)
  }
  total.list = list()
  n = length(methods.list.all)

  N = length(Y)
  B = length(levels(block))
  block.levels = levels(block)
  ## calculate numbers of observations in each block
  nb = rep(NA, B)
  for(i in 1:B){
    nb[i] = sum(block == block.levels[i])
  }


  for (l in 1:n){
    method.list.all = methods.list.all[[l]]

    Tlist = list()
    for(i in 1:B){
      Zb = Z[block == block.levels[i]]
      Yb = Y[block == block.levels[i]]
      Ti = matrix(nrow = 2, ncol = nb[i] + 1)
      Ti[1,] = 0:nb[i]
      for(ii in 0:nb[i]){
        if(length(method.list.all)==1){
          method.list = method.list.all[[1]]
        }else{
          method.list = method.list.all[[i]]
        }
        Ti[2, ii+1] = min_stat(Zb, Yb, nb[i]-ii, c, method.list = method.list)
      }
      Tlist[[i]] = Ti
    }
    total.list[[l]] <- Tlist
  }
  return(total.list)
}
