#' p_val_comb
#'
#' Calculate p-value for SCRE using combining multiple methods
#'
#' @export





###### p-value for testing tau_{(k)} \leq c, using combined statistic ######


pval_comb <- function(Z, Y, k, c,
                      block,
                      methods.list.all,
                      weight = NULL,
                      stat.null = NULL,
                      null.max = 10^5,
                      Z.perm.all = NULL,
                      statistic = TRUE,
                      switch = TRUE,
                      opt.method = "ILP_gurobi"){

  if(opt.method == "ILP_gurobi"){
    exact = TRUE
  } else {
    exact = FALSE
  }

  N = length(Y)
  p = N - k

  if(!is.factor(block)){
    block = as.factor(block)
  }

  
  if(switch){
    B = length(levels(block))
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
 
  ms_list = mu_sigma_list(Z = Z, block = block,
                          methods.list.all = methods.list.all)

  if(is.null(stat.null)){
    stat.null = com_null_dist_block(Z = Z, block = block,
                                    methods.list.all = methods.list.all,
                                    null.max = null.max,
                                    Z.perm.all = NULL,
                                    mu_sigma_list = ms_list)
  }

  coeflists = comb_matrix_block(Z = Z, Y = Y,
                                block = block, c = c,
                                methods.list.all = methods.list.all)

  stat.min = Gurobi_sol_com(Z = Z, block = block,
                            weight = weight,
                            coeflists = coeflists,
                            p = p,
                            mu_sigma_list = ms_list,
                            exact = exact)$obj
  pval = mean(stat.null >= stat.min)
  if(statistic == TRUE) {
    result = c(pval, stat.min)
    names(result) = c("p.value", "test.stat")
    return(result)
  } else return(pval)
}


