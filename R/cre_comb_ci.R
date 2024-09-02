### Last Update: 07/03/2024
### Codes for CMRSS on CRE
### Main Changes: (1) p-value is constructed by comparing minimum tail probabilities
###               (2)
###               (3) function for simultaneous inference is added(Zhe Chen and Xinran Li's method is also applied)

#library(CMRSS)

## function calculating one p-value for individual methods(including polynomial scores)

pval_cre <- function(Z, Y, k, c,
                     method.list, score = NULL, stat.null = NULL,
                     nperm = 10^3, Z.perm = NULL, ind.sort.treat = NULL){
  n = length(Z)
  m = sum(Z)

  # get score if score is null #
  if(is.null(score)){
    score = CMRSS::rank_score( n, method.list )
  }

  # emp null dist #
  if(is.null(stat.null)){
    stat.null = CMRSS::null_dist(n, m, score = score, nperm = nperm, Z.perm = Z.perm, method.list = method.list )
  }

  # min stat value under H_{k,c} #
  stat.min = min_stat(Z, Y, k, c, score = score, ind.sort.treat = ind.sort.treat, method.list = method.list )

  # p-value #
  pval = mean( stat.null >= stat.min )

  return(pval)
}

## function generating null distribution for minimum among tail probabilities

comb_null_dist_cre = function(n, m, methods.list, Z.perm = NULL,nperm = 10^4){

  H = length(methods.list)
  if(is.null(Z.perm)){
  Z.perm = RIQITE::assign_CRE(n, m, nperm)
  }
  nperm = ncol(Z.perm)
  stat.dist = matrix(NA, nrow = H, ncol = nperm)
  prob.dist = matrix(NA, nrow = H, ncol = nperm)

  for(i in 1 : H){
    method.list = methods.list[[i]]
    stat.dist[i,] = CMRSS::null_dist(n, m, method.list, Z.perm = Z.perm)
  }
  for(j in 1 : H){
  #  prob.dist[j,] = (nperm + 1 - rank(stat.dist[j,])) / nperm
    prob.dist[j,] = rank(stat.dist[j,]) / nperm
  }
#  for(j in 1 : H){
#    tmp.dist = stat.dist[j,]
#    for(iter in 1 : nperm){
#      prob.dist[j,iter] = mean(tmp.dist >= stat.dist[j,iter])
#    }
#  }

  result = apply(prob.dist, 2, min)
  return(result)
}

## function for calculating valid minimum p-value

min_p_val_cre = function(Z, Y, k, c, stat.null = NULL,
                         methods.list, nperm){

  H = length(methods.list)
  pval.vec = rep(NA, H)
  n = length(Z)
  m = sum(Z)


  Z.perm = RIQITE::assign_CRE(n, m, nperm)

  if(is.null(stat.null)){
    null.dist = comb_null_dist_cre(n = n, m = m,
                                 methods.list = methods.list, Z.perm = Z.perm, nperm = nperm)
  }else{
    null.dist = stat.null
  }

  for (i in 1 : H){
    method.list = methods.list[[i]]
    pval.vec[i] = pval_cre(Z = Z, Y = Y, k = k, c = c, Z.perm = Z.perm, method.list = method.list)
  }

  min.pval = min(pval.vec)
  result = mean(null.dist <= min.pval)
  return(result)
}

# Simultaneous Inference for multiple quantiles on CRE

## helper function
## maybe in the form of maximum test statistic will be much efficient, but not sure about their equivalence

com_conf_quant_larger_trt <- function( Z, Y, methods.list = NULL,
                                       stat.null = NULL,
                                       nperm = 10^4,
                                       k.vec = NULL,
                                       Z.perm = NULL,
                                       alpha = 0.05,
                                       tol = 10^(-3),
                                       ind.sort.treat = NULL ){
  n = length(Z)
  m = sum(Z)

  # get score if score is null #
  #  if(is.null(score)){
  #    score = rank_score( n, method.list )
  #  }

  # emp null dist #
  if(is.null(stat.null)){
    stat.null = comb_null_dist_cre(n, m, nperm = nperm, methods.list = methods.list )
    nperm = length(stat.null)
  }else{
    nperm = length(stat.null)
  }

  # find threshold such that p-value <= alpha <===> test stat > threshold #
  thres = sort(stat.null, decreasing = FALSE)[ floor(nperm * alpha) + 1 ]

  # range of c #
  Y1.max = max(Y[Z==1])
  Y1.min = min(Y[Z==1])
  Y0.max = max(Y[Z==0])
  Y0.min = min(Y[Z==0])

  c.max = Y1.max - Y0.min + tol
  c.min = Y1.min - Y0.max - tol


  # sort the treated units
  if(is.null(ind.sort.treat)){
    ind.sort.treat = sort_treat(Y, Z)
  }

  if( is.null(k.vec)) {
    # conf interval for all quantiles tau_{(k)} #
    c.limit = rep(NA, n)

    for(k in n:(n-m)){
      # define the target fun #
      # f > 0 <==> p-value <= alpha
      # f decreases in c, p value increases in c
      # we want to find minimum c which is close as possible to f <= 0; for 1-alpha CI
      if(k < n){
        c.max = c.limit[k+1]
      }
      f <- function(c){
        pval = min_p_val_cre(Z, Y, k, c, stat.null = stat.null, methods.list, nperm) # not a min_stat, but a min pval.
        return(thres - pval)
      }
      # check whether f(-Inf) = f(c.min) > 0 #
      if( f(c.min) <= 0 ){  ##
        c.sol = -Inf
      }
      if( f(c.max) > 0){  ##
        c.sol = c.max
      }
      if( f(c.min) > 0 & f(c.max) <= 0 ){   # need to change this in reverse (c.min, c.max)
        c.sol = uniroot(f, interval = c(c.min, c.max), extendInt = "downX", tol = tol)$root
        # find the min c st p-value > alpha <==> f <= 0 #
        c.sol = round(c.sol, digits = -log10(tol))
        if( f(c.sol) <= 0 ){
          while( f(c.sol) <= 0 ){
            c.sol = c.sol - tol
          }
          c.sol = c.sol + tol
        }else{
          while(f(c.sol) > 0 ){
            c.sol = c.sol + tol
          }
        }
      }
      c.limit[k] = c.sol
    }

    if( n-m > 1){
      c.limit[1:(n-m-1)] = c.limit[n-m]
    }

    #c.limit[c.limit > (Y1.max - Y0.min) + tol/2] = Inf
  } else {

    k.vec.sort = sort(k.vec, decreasing = FALSE)
    j.max = length(k.vec.sort)
    j.min = max( sum(k.vec <= (n-m)), 1)

    c.limit = rep(NA, j.max)

    for(j in j.max:j.min){
      k = k.vec.sort[j]
      if(j < j.max){
        c.max = c.limit[j+1]
      }
      f = function(c){
        pval = min_p_val_cre(Z, Y, k, c, stat.null, methods.list, nperm)
        return(thres - pval)
      }
      if( f(c.min) <= 0){
        c.sol = -Inf
      }
      if( f(c.max) > 0){
        c.sol = c.max
      }
      if(f(c.min) > 0 & f(c.max) <= 0){
        c.sol = uniroot(f, interval = c(c.min, c.max), extendInt = "downX", tol = tol)$root
        c.sol = round(c.sol, digits = -log10(tol))
        if( f(c.sol) <= 0){
          while( f(c.sol) <= 0){
            c.sol = c.sol - tol
          }
          csol = c.sol + tol
        }else{
          while(f(c.sol) > 0 ){
            c.sol = c.sol + tol
          }
        }
      }
      c.limit[j] = c.sol
    }
    if(j.min > 1){
      c.limit[1:(j.min-1)] = c.limit[j.min]
    }
#    c.limit[c.limit > (Y1.max - Y0.min) + 10^{-log10(tol)} / 2] = Inf
  }

  return( c.limit )
}

## main function
## generates lower bound of each confidence interval of c's on every k's. (so for all individual treatment effects)

com_conf_quant_larger_cre = function( Z, Y, methods.list = NULL,
                                      stat.null = NULL, # by now, this is just for the treated ones(will add option for scores for control units)
                                      nperm = 10^4,
                                      set = "prev",
                                      Z.perm = NULL,
                                      k.vec = NULL,
                                      alpha = 0.05,
                                      tol = 10^(-3),
                                      ind.sort.treat = NULL ){

  H = length(methods.list)
  n = length(Z)
  # set = "prev": without Chen and Li's adjustment
  if(set == "prev"){
    ci.prev = com_conf_quant_larger_trt(Z, Y, methods.list, stat.null,
                                        nperm, k.vec,
                                        Z.perm, alpha, tol, ind.sort.treat
    )
    return(ci.prev)
  }

  if(set == "trt"){
    ci.treat = com_conf_quant_larger_trt(Z, Y, methods.list, stat.null,
                                         nperm, k.vec = (n - sum(Z)+ 1) : n,
                                         Z.perm, alpha, tol, ind.sort.treat
    )
    return(ci.treat)
  }

  if(set == "control"){
    stat.null.control = NULL
    Z.perm.control = NULL
    Z = 1 - Z
    Y = -Y
    ci.control = com_conf_quant_larger_trt(Z, Y,
                                           methods.list = methods.list,
                                           stat.null = stat.null.control,
                                           nperm = nperm,
                                           k.vec = (n - sum(Z) + 1) : n,
                                           Z.perm = Z.perm.control,
                                           alpha = alpha, tol = tol
    )
    return(ci.control)
  }

  # set = "all": with Chen and Li's adjustment
  if(set == "all"){
    # trted ones
    ci.treat = com_conf_quant_larger_trt(Z, Y, methods.list, stat.null,
                                         nperm, k.vec = (n - sum(Z) + 1):n,
                                         Z.perm, alpha, tol, ind.sort.treat)

    # controlled ones
    stat.null.control = NULL
    Z.perm.control = NULL
    Z = 1 - Z
    Y = -Y
    ci.control = com_conf_quant_larger_trt(Z, Y,
                                           methods.list = methods.list,
                                           stat.null = stat.null.control,
                                           nperm = nperm,
                                           k.vec = (n - sum(Z) + 1) : n,
                                           Z.perm = Z.perm.control,
                                           alpha = alpha, tol = tol
    )
    ci.all = sort( c(ci.treat, ci.control))
    return(ci.all)
  }
}


######################################################################
# double check function works / comparing to original function in RIQITE
#######################################################################

#library(mvtnorm)
#n = 120
#m = 0.5 * n
#Z_block = factor(rep(1, n))
#nperm = 10^3
#alpha = 0.1
#k = n #ceiling(n * 0.9)
#p = n - k
#c = 0

#s.vec = floor(seq(from = 2, to = m, len = 6))
#H = length(s.vec)

#methods.list.com = list()
#for(i in 1 : H){
#  methods.list.com[[i]] = list()
#}
#for(i in 1 : H){
#  methods.list.com[[i]]           = list(name = "Polynomial",
#                                         r = s.vec[i],
#                                         std = TRUE,
#                                         scale = TRUE)
#}
#method.list = list(name = "Stephenson", s = 32)

#w = 0.5
#tau_0 = 0.5
#mu = c(0, tau_0)

#sigma = matrix(c(1, 0.5 * w, 0.5 * w, w^2), nrow = 2)

#tmp = rmvnorm(n = n, mean = mu, sigma = sigma)

#Y0 = tmp[,1]
#tau = tmp[,2]
#Y1 = Y0 + tau

#Z = assign_CRE(n, m, 1)
#Y = Z * Y1 + (1-Z) * Y0


# check individual (combined) p-values is printed

#pval_cre(Z, Y, k, c, method.list = methods.list.com[[3]])

#min_p_val_cre(Z, Y, k, c,
#              methods.list = methods.list.com, nperm = 10^4)

# simultaneous conf. interval on CRE

#com_conf_quant_larger_cre(Z, Y, methods.list = methods.list.com, nperm = 10^3, set = "trt")
#com_conf_quant_larger_cre(Z, Y, methods.list = methods.list.com, nperm = 10^3, set = "all") # applying Chen and Li's adjustment

## compare the result with RIQITE
#library(RIQITE)
#pval_quantile(Z, Y, k, c, method.list = methods.list.com[[3]], nperm = nperm)

