### Last Update: 08/16/2024
### function for simultaneous inference on the quantile of individual treatment effects, on SERE
### generates lower bound of each confidence interval of c's on every k's.


#' calculating null distribution on cmrss

com_null_dist_block <- function(Z, block,
                                methods.list.all,
                                null.max = 10^5,
                                Z.perm.all = NULL,
                                ms.list){
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
      temp[j,] = temp[j,] + CMRSS::null_dist(nb[i], mb[i],
                                             method.list = method.list,
                                             Z.perm = Z.perm.i)
    }
  }

  for(i in 1 : K){
    temp[i,] = (temp[i,] - ms.list$mean[i]) / ms.list$sigma[i]
  }
  stat.null = apply(temp, 2, max)
  return(stat.null)
}

#' helper function
com_block_conf_quant_larger_trt = function(Z, Y,
                                           block,
                                           quantiles = NULL,
                                           methods.list.all = NULL,
                                           opt.method = "ILP_gurobi",
                                           stat.null = NULL, null.max = 10^4,
                                           tol = 0.01,
                                           alpha = 0.1){



  if(!is.factor(block)){
    block = as.factor(block)
  }
  levels(block) = 1:length(levels(block))
  n = length(Z)
  m = sum(Z)
  B = length(levels(block))
  ms.list = mu_sigma_list(Z, block, methods.list.all)

  if( all(is.nan(ms.list$mean)) ) {
    stop("Check choice of test parameter is larger than the block size")
  }

  if(is.null(stat.null)){
    stat.null = com_null_dist_block(Z, block, methods.list.all,
                                    null.max, ms.list = ms.list)
  }
  thres = sort(stat.null, decreasing = TRUE)[floor(null.max * alpha) + 1]

  Y1.max = max(Y[Z == 1])
  Y0.min = min(Y[Z == 0])

  c.max = max(Y[Z==1]) - min(Y[Z == 0]) + tol
  c.min = min(Y[Z==1]) - max(Y[Z == 0]) - tol

  #### should we do the ind.sort.treat here? => I think this is already added in min_stat from CMRSS

  if( is.null(quantiles) ){

    quantiles = rep(NA, n)

    for(k in n : (n-m)){

      if(k < n){
        c.max = quantiles[k+1]
      }

      f = function(c){
        stat.min = pval_comb(Z, Y, k, c,
                             block,
                             methods.list.all,
                             null.max = null.max,
                             statistic = TRUE,
                             switch = FALSE)[2]
        return(stat.min - thres)
      }
      if( f(c.min) <= 0){
        c.sol = -Inf
      }
      if( f(c.max) > 0){
        c.sol = c.max
      }
      if( f(c.min) > 0 & f(c.max) <= 0 ){
        c.sol = uniroot(f, interval = c(c.min, c.max), extendInt = "downX", tol = tol)$root
        c.sol = round(c.sol, digits = -log10(tol))

        if( f(c.sol) <= 0 ){
          while( f(c.sol) <= 0){
            c.sol = c.sol - tol
          }
          c.sol = c.sol + tol
        }else{
          while(f(c.sol) > 0){
            c.sol = c.sol + tol
          }
        }
      }
      quantiles[k] = c.sol
    }

    if( n-m > 1){
      quantiles[1 : (n-m-1)] = quantiles[n-m]
    }

    quantiles[quantiles > (Y1.max - Y0.min) + tol/2] = Inf

  }else{

    k.vec.sort = sort(quantiles, decreasing = FALSE)
    j.max = length(k.vec.sort)
    j.min = max( sum(k.vec <= (n-m)), 1)
    quantiles = rep(NA, j.max)

    for(j in j.max:j.min){
      k = k.vec.sort[j]
      if(j < j.max){
        c.max = quantiles[j+1]
      }
      f = function(c){
        stat.min = pval_comb(Z, Y, k, c,
                             block,
                             methods.list.all,
                             null.max = null.max,
                             statistic = TRUE,
                             switch = FALSE)[2]
        return(stat.min - thres)
      }
      if( f(c.min) <= 0){
        c.sol = -Inf
      }
      if( f(c.max) > 0){
        c.sol = c.max
      }
      if( f(c.min > 0) & f(c.max) <= 0){
        c.sol = uniroot(f, interval = c(c.min, c.max), extendInt = "downX", tol = tol)$root
        c.sol = round(c.sol, digits = -log10(tol))
        if( f(c.sol) <= 0){
          while( f(c.sol) <= 0){
            c.sol = c.sol - tol
          }
          c.sol = c.sol + tol
        }else{
          while(f(c.sol) > 0){
            c.sol = c.sol + tol
          }
        }
      }
      quantiles[j] = c.sol
    }
    if(j.min > 1){
      quantiles[1 : (j.min - 1)] = quantiles[j.min]
    }

  }
  return(quantiles)
}




#' Simultaneous confidence interval using CMRSS on SERE
#'
#' @export

com_block_conf_quant_larger = function(Z, Y,
                                       block,
                                       quantiles = NULL,
                                       set = "all",
                                       methods.list.all = NULL,
                                       opt.method = "ILP_gurobi",
                                       stat.null = NULL, null.max = 10^4, alpha = 0.1){
  if(set == "control"){
    Z = 1 - Z
    Y = -Y
  }

  if(set == "all"){
    Z1 = 1 - Z
    Y1 = -Y
  }

  ci.1 = com_block_conf_quant_larger_trt(Z, Y, block, quantiles, methods.list.all, opt.method,
                                         stat.null, null.max, alpha)
  if(set == "trt" | set == "control"){
    return(as.numeric(ci.1))
  }else{
    ci.2 = com_block_conf_quant_larger_trt(Z1, Y1, quantiles, methods.list.all, opt.method,
                                           stat.null, null.max, alpha )
    ci.2 = as.numeric(ci.2)
    ci.all = sort(c(ci.1, ci.2))
    return(ci.all)
  }
}


# double check function works / comparing to original function in QIoT

#source("sere_data_generate.R")
#library(QIoT)

#ci_quantile_scre(Z, Y, Z_block,
#                   method.list.all = method.list.all, opt.method = "ILP_gurobi", null.max = 10^3)

#com_block_conf_quant_larger_trt(Z, Y, block,
#                                methods.list.all = methods.list.all,
#                                null.max = 10^3, alpha = 0.1)








