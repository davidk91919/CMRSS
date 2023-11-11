#' com_null_dist_block
#'
#' @param Z n-dimension treatment assignment vector, for every units
#' @param block Length n vector specifies stratum status for every units
#' @param method.list.all list of every stratified rank sum statistics, where each element(list)s is method for every stratum.List of stratified rank sum statistic. The structure should be list of lists of lists.
#' @param null.max Total size of randomization distribution.
#' @param Z.perm.all User-specified assignment matrix, for every element of null distribution.
#' @param mu_sigma_list mu_sigma_list List of population mean / variance of every stratified rank sum statistics. This will usually be a output from \emph{mu_sigma_list} function.
#' @noRd


com_null_dist_block <- function(Z, block,
                                methods.list.all,
                                null.max = 10^5,
                                Z.perm.all = NULL,
                                mu_sigma_list){
  if(!is.factor(block)){
    block = as.factor(block)
  }

  K = length(methods.list.all)

  N = length(Z)
  B = length(levels(block))
  b_num = as.vector(table(block))
  block.levels = levels(block)
  temp = matrix(0, nrow = K, ncol = null.max)

  nb = rep(NA, B)
  mb = rep(NA, B)
  for(i in 1:B) {
    nb[i] = sum(block == block.levels[i])
    mb[i] = sum( Z[block == block.levels[[i]]] )
  }

  for(i in 1:B){
    if(is.null(Z.perm.all)){
      Z.perm.i = matrix(nrow = nb[i], ncol = null.max)
      for(iter in 1:null.max){
        Z.perm.i[,iter] = sample( c(rep(1, mb[i]), rep(0, nb[i] - mb[i])) )
      }
    }else{
      Z.perm.i[,iter] = Z.perm.all[block == block.levels[i], ]
    }
    for(j in 1 : K){
      method.list.all = methods.list.all[[j]]

      if(length(method.list.all) == 1){
        method.list = method.list.all[[1]]
      }else{
        method.list = method.list.all[[i]]
      }
      temp[j,] = temp[j,] + null_dist(nb[i], mb[i],
                                      method.list = method.list,
                                      Z.perm = Z.perm.i)
    }
  }

  for(i in 1 : K){
    temp[i,] = (temp[i,] - mu_sigma_list$mean[i]) / mu_sigma_list$sigma[i]
  }
  stat.null = apply(temp, 2, max)
  return(stat.null)
}



