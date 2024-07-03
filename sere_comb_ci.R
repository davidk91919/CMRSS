
# function for simultaneous inference on the quantile of individual treatment effects

library(CMRSS)

## helper function

com_block_conf_quant_larger_trt = function(Z, Y,
                                       block,
                                       quantiles = NULL,
                                       methods.list.all = NULL,
                                       opt.method = "ILP_gurobi",
                                       stat.null = NULL, null.max = 10^4, alpha = 0.1){
  if(!is.factor(block)){
    block = as.factor(block)
  }
  levels(block) = 1:length(levels(block))
  n = length(Z)
  mn = sum(Z)
  B = length(levels(block))
  
  if(switch){
    nb = rep(NA, B)
    zb = rep(NA, B)
    for(i in 1:B){
      nb[i] = sum(block == levels(block)[i])
      zb[i] = sum(Z[block == levels(block)[i]])
    }
    for(i in 1:B){
      if(zb[i]<nb[i]/2){
        Z[block == levels(block)[i]] = 1 - Z[block == levels(block)[i]]
        Y[block == levels(block)[i]] = -Y[block == levels(block)[i]]
      }
    }
  }
  
  ms_list = mu_sigma_list(Z, block, methods.list.all = methods.list.all)
  
  if(is.null(stat.null)){
    stat.null = com_null_dist_block(Z, block,
                                    methods.lists.all = methods.lists.all,
                                    null.max = null.max, mu_sigma_list = ms_list)
  }
  
  thres = sort(stat.null, decreasing = TRUE)[ floor(length(stat.null) * alpha) + 1]
  cup = max(Y[Z == 1]) - min(Y[Z == 0]) + 0.01
  cdown = min(Y[Z == 1]) - max(Y[Z == 0]) - 0.01
  cmid = (cup + cdown) / 2
  
  if(is.null(quantiles)){
    quantiles = 1 : n   # different with RIQITE; they used NA for the initial vector for quantiles.
  }
  
  
  
  quantiles = sort(quantiles)
  c_conf1 = rep(NA, length(quantiles))
  
  for (k in length(quantiles):1){
    coeflists = comb_matrix_block(Z, Y, block, cdown, methods.list.all = methods.list.all)
    up = Gurobi_sol_com(Z, block, p = n - k, coeflists = coeflists, mu_sigma_list = ms_list, exact = TRUE)$obj
    if(up > thres & quantiles[k] != n){
      c_conf1[k] = cmid
      next
    }
    
    if(up < thres){
      c_conf1[1:k] = -Inf
      break
    }
    else
    {
      repeat{
        coeflists = comb_matrix_block(Z, Y, block, cmid, methods.list.all = methods.list.all)
        mid = Gurobi_sol_com(Z, block, p = n - k, coeflists = coeflists, mu_sigma_list = ms_list)$obj
        if(mid > thres){
          cdown = cmid
        }else{
          cup = cmid
        }
        
        if(abs(cup-cdown) < 1e-6){
          c_conf1[k] = cmid
          break
        }
      }
    }
  }
  return(c_conf1)
}

## main function 

com_block_conf_quant_larger_trt = function(Z, Y,
                                           block,
                                           quantiles = NULL,
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


