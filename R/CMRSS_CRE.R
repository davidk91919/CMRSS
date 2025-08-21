
#' sorting index among treated units, according to ranks of the outcome
#'
#' output: indices of treated units with increasing outcomes

sort_treat <- function(Y, Z){
  r = rank(Y, ties.method = "first")
  ind.sort = sort.int(r, index.return = TRUE)$ix
  ind.sort.treat = ind.sort[Z[ind.sort] == 1]
  return(ind.sort.treat)
}



#' Generate matrix of complete randomization assignments
#' output: a matrix with nperm columns, each column gives one assignment
assign_CRE <- function(n, m, nperm){
  if(is.finite(nperm)){
    Z.perm = matrix(0, nrow = n, ncol = nperm)
    for(iter in 1:nperm){
      Z.perm[sample( c(1:n), m, replace = FALSE ), iter] = 1
    }
  }

  if(is.infinite(nperm)){
    comb.all = combn(n, m)
    nperm = ncol(comb.all)
    Z.perm = matrix(0, nrow = n, ncol = nperm)
    for(iter in 1:nperm){
      Z.perm[comb.all[, iter], iter] = 1
    }
  }

  return(Z.perm)
}



#' Calculating score for given n
#' included polynomial rank score with Puri(1965)'s normalization (denoted by std)
#' output: rank scores of length n
rank_score <- function(n, method.list = list(name = "Polynomial", r, std = TRUE, scale = FALSE) ){
  if(method.list$name == "Polynomial"){
    r = method.list$r
    if(method.list$std == TRUE) {
      #      score = (c(1:n) / n)^(r-1)
      score = (c(1:n) / (n + 1) )^(r-1)  # changed
    } else {
      score = (c(1:n))^(r-1)
    }
    if(method.list$scale == TRUE){
      score = scale(score) # normalized score; corresponds to 11p of the draft
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
  return(score)
}



#' generating randomization null distribution of a single rank sum statistic
#' output: null dist of a single rank score statistics under the CRE
null_dist <- function(n, m, method.list = NULL, score = NULL,
                      nperm = 10^5, Z.perm = NULL){
  if(is.null(score)){
    score = rank_score( n, method.list )
  }

  if(is.null(Z.perm)){
    Z.perm = assign_CRE(n, m, nperm)
    nperm = ncol(Z.perm)
  }

  stat.null = rep(NA, ncol(Z.perm))
  for(iter in 1:ncol(Z.perm)){
    stat.null[iter] = sum( score[ Z.perm[, iter] == 1 ] )
  }

  return(stat.null)
}

#' null distribution of multiple rank sum statistics
#' output: matrix of null distributions from multiple rank statistics
null_dist_multiple <- function(n, m, methods.list = NULL, nperm = 10^5, Z.perm = NULL){
  if(is.null(Z.perm)){
    Z.perm = assign_CRE(n, m, nperm)
    nperm = ncol(Z.perm)
  }

  H = length(methods.list)
  stat.null.mult = matrix(NA, nrow = H, ncol = nperm)
  for(h in 1:H){
    stat.null.mult[h, ] = null_dist(n, m, method.list = methods.list[[h]], score = NULL, nperm = nperm, Z.perm = Z.perm)
  }

  return(stat.null.mult)
}



#' calculating individual test statistic corresponds to Caughey et al.(2023)
#' output: minimum value of the rank sum statistic assuming n-k treated units can have infinite effects
min_stat <- function(Z, Y, k, c, method.list = NULL,
                     score = NULL,
                     ind.sort.treat = NULL){
  n = length(Y)
  m = sum(Z)

  # if(is.null(method.list)){   #added # I delete, we can have scale in the method.list
  #   score= scale(score)
  # }

  if(is.null(score)){
    score = rank_score(n, method.list)
  }

  # sort the treated units
  if(is.null(ind.sort.treat)){
    ind.sort.treat = sort_treat(Y, Z)
  }

  # get xi vector
  xi = rep(c, n)
  if(k < n){
    xi[ ind.sort.treat[ ( m + 1 - min(m, n-k) ):m ] ] = Inf
  }

  stat.min = sum( score[ rank( Y - Z * xi, ties.method = "first" )[Z==1] ] )

  return(stat.min)
}



#' calculating p-value for single method
#' output: the p-value assuming at most n-k treated units having effects > c, using a single rank sum statistic
pval_cre <- function(Z, Y, k, c,
                     method.list, score = NULL, stat.null = NULL,
                     nperm = 10^3, Z.perm = NULL, ind.sort.treat = NULL){
  n = length(Z)
  m = sum(Z)

  # get score if score is null #
  if(is.null(score)){
    score = rank_score( n, method.list )
  }

  # emp null dist #
  if(is.null(stat.null)){
    if(is.null(Z.perm)){
      Z.perm = assign_CRE(n, m, nperm)
      nperm = ncol(Z.perm)
    }
    stat.null = null_dist(n, m, score = score, nperm = nperm, Z.perm = Z.perm, method.list = method.list )
  }

  # min stat value under H_{k,c} #
  stat.min = min_stat(Z, Y, k, c, score = score, ind.sort.treat = ind.sort.treat, method.list = method.list )

  # p-value #
  pval = mean( stat.null >= stat.min )

  return(pval)
}



#' null distribution for minimum p-values
#' output: the distribution of the minimum p-value, particularly Monte Carlo samples from \overline{F} in Theorem 1 of the paper
comb_null_dist_cre = function(n, m, methods.list, Z.perm = NULL, nperm = 10^4, stat.null.mult = NULL){

  H = length(methods.list)

  if(is.null(stat.null.mult)){
    if(is.null(Z.perm)){
      Z.perm = assign_CRE(n, m, nperm)
    }
    stat.null.mult = null_dist_multiple(n, m, methods.list, nperm, Z.perm)
  }

  tail.prob = matrix(NA, nrow = H, ncol = nperm)
  for(j in 1 : H){
    tail.prob[j, ] = 1 - (rank(stat.null.mult[j,], ties.method = "min") - 1) / nperm
  }

  result = apply(tail.prob, 2, min)
  return(result)
}



#' calculate the minimum p value from multiple rank sum statistics
#' output: minimum of p values calibrated by the distribution of minimum of tail probabilities
min_p_multiple_rank_sum <- function(Z, Y, k, c, methods.list, Z.perm = NULL, nperm, stat.null.mult = NULL){
  
  n = length(Z)
  m = sum(Z)
  
  if(is.null(stat.null.mult)){
    if(is.null(Z.perm)){
      Z.perm = assign_CRE(n, m, nperm)
    }
    stat.null.mult = null_dist_multiple(n, m, methods.list, nperm, Z.perm)
  }

  H = length(methods.list)
  pval.vec = rep(NA, H)
  for (i in 1 : H){
    method.list = methods.list[[i]]
    pval.vec[i] = pval_cre(Z = Z, Y = Y, k = k, c = c, Z.perm = Z.perm, method.list = method.list, stat.null = stat.null.mult[i, ])
  }

  min.pval = min(pval.vec)
  return(min.pval)
}


#' Function for calculating valid minimum p-value in CRE
#' Calculate a valid p-value, based on multiple random sum statistics,
#' for testing the null that at most n-k treated units having effects > c,
#' null hypothesis H0: \eqn{\tau_{(k)} \leq c}, or \eqn{H0: \tau_{(k)} \geq c}, or
#' \eqn{H0: \tau_{(k)} = c}, where \eqn{\tau_{(k)}} denotes individual
#' treatment effect at rank k.
#'
#' @param Z An \eqn{n} dimensional treatment assignment vector.
#' @param Y An \eqn{n} dimensional outcome vector.
#' @param k An integer between 1 and n specifying which quantile of
#'   individual effect is of interest.
#' @param c A numerical object specifying the threshold for the null
#'   hypothesis.
#' @param methods.list A list of lists specifies the choice of the multiple
#'   rank sum test statistics. For example, list(name = "Wilcoxon") means the Wilcoxon
#'   rank sum statistic, list(name = "Stephenson", s = 10) means
#'   the Stephenson rank sum statistic with parameter s = 10, and
#'   list(name = "Polynomial", r = 2, std = TRUE, scale = FALSE) means using polynomial
#'   rank score with standarizing statistics without scaling.
#' @param Z.perm A \eqn{n \times nperm} matrix that specifies the
#'   permutated assignments for approximating the null distribution of
#'   the test statistic.
#' @param nperm A positive integer representing the number of
#'   permutations for approximating the randomization distribution of
#'   the rank sum statistic.
#' @param stat.null.mult A matrix whose empiricial distribution
#'   approxmiates the randomization distribution of multiple rank statistics.
#' @param null.dist.comb
#'
#' @return Combined p-value for testing the specified null hypothesis of
#'   interest in CRE.
#' @export

comb_p_val_cre = function(Z, Y, k, c, methods.list,
                          Z.perm = NULL, nperm,
                          stat.null.mult = NULL
                          ){
  H = length(methods.list)
  pval.vec = rep(NA, H)
  n = length(Z)
  m = sum(Z)

  if(is.null(stat.null.mult)){
    if(is.null(Z.perm)){
      Z.perm = assign_CRE(n, m, nperm)
    }
    stat.null.mult = null_dist_multiple(n, m, methods.list, nperm, Z.perm)
  }

  null.dist.comb = comb_null_dist_cre(n, m, methods.list, Z.perm, nperm, stat.null.mult)

  min.pval = min_p_multiple_rank_sum(Z, Y, k, c, methods.list, Z.perm, nperm, stat.null.mult)
  result = mean(null.dist.comb <= min.pval) # XL: I changed
  return(result)
}


#' Helper function for Simultaneous Inference for multiple quantiles on CRE
#' using combined p-value.
#'
#' Output: lower limits of prediction intervals for (prespecified) quantiles across the treated units.
com_conf_quant_larger_trt <- function( Z, Y, methods.list = NULL,
                                       nperm = 10^4,
                                       k.vec = NULL,
                                       Z.perm = NULL,
                                       alpha = 0.05,
                                       tol = 10^(-3),
                                       ind.sort.treat = NULL ){
  n = length(Z)
  m = sum(Z)

  # emp null dist #
  # XL: I add the following
  if(is.null(Z.perm)){
    Z.perm = assign_CRE(n, m, nperm)
  }else{
    nperm = ncol(Z.perm)
  }

  stat.null.mult = null_dist_multiple(n, m, methods.list, nperm, Z.perm)
  null.dist.comb = comb_null_dist_cre(n, m, methods.list, Z.perm, nperm, stat.null.mult)

  # find threshold such that valid p-value <= alpha <===> min_pval < something <===>  -min_pval > -something #
  thres.minus.min.pval = sort(-1*null.dist.comb, decreasing = TRUE)[ floor(nperm * alpha) + 1 ]
  thres = alpha

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
      ## XL: an alternative way that using minus min pval as the test statistic
      f <- function(c){
        stat.min = -1 * min_p_multiple_rank_sum(Z, Y, k, c, methods.list, Z.perm, nperm, stat.null.mult)
        return(stat.min - thres.minus.min.pval)
      }

      # check whether f(-Inf) = f(c.min) > 0 #
      if( f(c.min) <= 0 ){
        c.sol = -Inf
      }
      if( f(c.max) > 0){
        c.sol = c.max
      }
      if( f(c.min) > 0 & f(c.max) <= 0 ){
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
  }else{

    k.vec.sort = sort(k.vec, decreasing = FALSE)
    j.max = length(k.vec.sort)
    j.min = max( sum(k.vec <= (n-m)), 1)

    c.limit = rep(NA, j.max)

    for(j in j.max:j.min){
      k = k.vec.sort[j]
      if(j < j.max){
        c.max = c.limit[j+1]
      }
      ## XL: an alternative way that using minus min pval as the test statistic
      f <- function(c){
        stat.min = -1 * min_p_multiple_rank_sum(Z, Y, k, c, methods.list, Z.perm, nperm, stat.null.mult)
        return(stat.min - thres.minus.min.pval)
      }

      # check whether f(-Inf) = f(c.min) > 0 #
      if( f(c.min) <= 0 ){
        c.sol = -Inf
      }
      if( f(c.max) > 0){
        c.sol = c.max
      }
      if( f(c.min) > 0 & f(c.max) <= 0 ){
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

      c.limit[j] = c.sol
    }
    if(j.min > 1){
      c.limit[1:(j.min-1)] = c.limit[j.min]
    }
  }

  return( c.limit )
}


#' Simultaneous lower bounds using multiple rank sum statistics in CRE
#'
#'
#' @param Z An \eqn{n} dimensional treatment assignment vector.
#' @param Y An \eqn{n} dimensional observed outcome vector.
#' @param methods.list A list of lists specifies the choice of the multiple
#'   rank sum test statistics. For example, list(name = "Wilcoxon") means the Wilcoxon
#'   rank sum statistic, list(name = "Stephenson", s = 10) means
#'   the Stephenson rank sum statistic with parameter s = 10, and
#'   list(name = "Polynomial", r = 2, std = TRUE, scale = FALSE) means using polynomial
#'   rank score with standarizing statistics without scaling.
#' @param nperm A positive integer representing the number of
#'   permutations for approximating the randomization distribution of
#'   the rank sum statistic.
#' @param set set of quantiles of interests. If "treat", generate prediction intervals
#'   for effect quantiles among treated units. If "control", generate prediction intervals
#'   for effect quantiles among control units. If "all", generate confidence intervals for
#'   all effect quantiles.
#' @param Z.perm A \eqn{n \times nperm} matrix that specifies the
#'   permutated assignments for approximating the null distribution of
#'   the test statistic.
#' @param alpha A numerical object, where 1-alpha indicates the
#'   confidence level.
#' @param tol A numerical object specifying the precision of the
#'   obtained confidence intervals. For example, if tol = 10^(-3),
#'   then the confidence limits are precise up to 3 digits.
#'
#' @return A vector specifying lower limits of prediction (confidence) intervals for
#' quantiles k = 1 ~ m (or n - m, or n).
com_conf_quant_larger_cre <- function( Z, Y, methods.list,
                                       nperm = 10^4,
                                       set = "treat",
                                       Z.perm = NULL,
                                       alpha = 0.05,
                                       tol = 10^(-3)
){

  H = length(methods.list)
  n = length(Z)
  if(is.null(Z.perm)){
    Z.perm = assign_CRE(n, m=sum(Z), nperm)
  }

  # prediction intervals for effect quantiles among treated units
  if(set == "treat"){
    ci.treat = com_conf_quant_larger_trt(Z, Y, methods.list,
                                         nperm,
                                         k.vec = NULL,
                                         Z.perm, alpha, tol, ind.sort.treat = NULL
    )
    ci.treat = ci.treat[(n - sum(Z) + 1):n]
    return(ci.treat)
  }

  # prediction intervals for effect quantiles among control units
  if(set == "control"){
    # stat.null.control = NULL
    if(is.null(Z.perm)){
      Z.perm.control = NULL
    }else{
      Z.perm.control = 1 - Z.perm
    }

    Z = 1 - Z
    Y = -Y
    ci.control = com_conf_quant_larger_trt(Z, Y,
                                           methods.list = methods.list,
                                           nperm = nperm,
                                           k.vec = NULL,
                                           Z.perm = Z.perm.control,
                                           alpha = alpha, tol = tol, ind.sort.treat = NULL)
    ci.control = ci.control[(n - sum(Z) + 1):n]
    return(ci.control)
  }

  # set = "all": with Chen and Li's adjustment
  if(set == "all"){
    # trted ones
    ci.treat = com_conf_quant_larger_trt(Z, Y, methods.list,
                                         nperm, k.vec = NULL,
                                         Z.perm, alpha = alpha/2, tol, ind.sort.treat=NULL)
    ci.treat = ci.treat[(n - sum(Z) + 1):n]

    # controlled ones
    # stat.null.control = NULL

    if(is.null(Z.perm)){
      Z.perm.control = NULL
    }else{
      Z.perm.control = 1 - Z.perm
    }
    Z = 1 - Z
    Y = -Y
    ci.control = com_conf_quant_larger_trt(Z, Y,
                                           methods.list = methods.list,
                                           nperm = nperm,
                                           k.vec = NULL,
                                           Z.perm = Z.perm.control,
                                           alpha = alpha/2, tol = tol, ind.sort.treat=NULL)
    ci.control = ci.control[(n - sum(Z) + 1):n]

    ci.all = sort( c(ci.treat, ci.control))
    return(ci.all)
  }
}
