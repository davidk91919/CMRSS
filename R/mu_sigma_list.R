#' mu_sigma_list
#'
#' @param Z n-dimensional treatment assignment vector, for every units.
#' @param block n-dimensional vector, which specifies stratum status of each units.
#' @param methods.list.all List of stratified rank sum statistic. The structure should be list of lists of lists.
#' @noRd



##############
# function for calculating population mean / variance for each statistics
##############

mu_sigma_list = function(Z, block, methods.list.all){
  if(!is.factor(block)){
    block = as.factor(block)
  }

  H = length(methods.list.all)
  mean = rep(0, H)
  sigma = rep(0, H)

  #  B = s = length(methods.list.all[[1]][[1]])
  B = s = length(methods.list.all[[1]])
  block.levels = levels(block)
  nb = rep(NA, B)
  for (i in 1 : B){
    nb[i] = sum(block == block.levels[i])
  }


  for(i in 1 : H){
    method.list.all = methods.list.all[[i]]
    tmp1 <- 0
    tmp2 <- 0
    for(j in 1 : s){
      method.list = method.list.all[[j]]
      ns = nb[j]
      ms = sum(Z[block == block.levels[j]])

      score = rank_score(ns, method.list)
      tmp1 = tmp1 + ms * 1/ns * sum(score)
      tmp2 = tmp2 + ms^2 * (1/ms - 1/ns) * (1 / (ns - 1)) *
        sum( (score - 1 / ns * sum(score))^2 )
      #       + (sum(score^2) - 1 / ns * sum(score)^2)      #same
    }
    mean[i] = tmp1
    sigma[i] = sqrt(tmp2)
  }
  return(list(mean = mean, sigma = sigma))
}
