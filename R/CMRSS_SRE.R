
#' function giving summary of the block
#' 
#' output is summary of the given SRE 
summary_block <- function(Z, block){
  # number of blocks
  if(!is.factor(block)){
    block = as.factor(block)
  }
  block.levels = levels(block)
  B = length(block.levels)
  
  # number of treated and control within each block
  nb = rep(NA, B)
  mb = rep(NA, B)
  units.block = list()
  for(i in 1 : B){
    units.block[[i]] = which(block == block.levels[i])
    nb[i] = length(units.block[[i]])
    mb[i] = sum(Z[units.block[[i]]])
  }
  mb_ctrl = nb - mb
  
  result = list(block=block, B = B, nb = nb, mb = mb, mb_ctrl = mb_ctrl, units.block = units.block, block.levels = block.levels )
  
  return(result)
}

#' function assigning treatment assignment for given summary of SRE
#' 
#' output n-length treatment assignment vector regarding SRE. 
assign_block <- function(block.sum, null.max = 10^4){
  block = block.sum$block
  B = block.sum$B
  nb = block.sum$nb
  mb = block.sum$mb
  mb_ctrl = block.sum$mb_ctrl
  units.block = block.sum$units.block
  block.levels = block.sum$block.levels
  
  n = sum(nb)
  
  Z.perm = matrix(NA, nrow = n, ncol = null.max)
  for(iter in 1 : null.max){
    tmp = rep(0, n)
    for(i in 1 : B){
      tmp[units.block[[i]]] = sample( c(rep(1, mb[i]), rep(0, mb_ctrl[i])) )
    }
    Z.perm[,iter] = tmp
  }
  
  return(Z.perm)
}

#' giving the weight under given scheme
#' 
#' output is a vector of giving weight for prespecified option of weights.
weight_scheme <- function(block.sum, weight.name = "asymp.opt"){
  B = block.sum$B
  nb = block.sum$nb
  mb = block.sum$mb
  mb_ctrl = block.sum$mb_ctrl
  
  weight = rep(NA, B)
  
  if(weight.name == "asymp.opt"){
    weight = mb
  }
  if(weight.name == "dis.free"){
    weight = (nb + 1) / mb_ctrl
  }
  
  return(weight)
}

#' score for each of the blocks under a given method.list.all
#' 
#' output is a list of score output for each blocks, when a method.list.all is given.
#' method.list.all is a list of lists
score_all_blocks <- function(nb, method.list.all){
  B = length(method.list.all)
  score.list.all = list()
  for(i in 1 : B){
    method.list = method.list.all[[i]]
    score = rank_score(nb[i], method.list)
    score.list.all[[i]] = score
  }
  return(score.list.all)
}


#' function to calculate single stratified rank sum statistic with weights
#' output is a scalar.

single_weight_strat_rank_sum_stat <- function(Z, Y, block, method.list.all, score.list.all = NULL, weight, block.sum = NULL){
  
  if(is.null(block.sum)){
    block.sum = summary_block(Z, block)
  }
  block = block.sum$block
  B = block.sum$B
  nb = block.sum$nb
  mb = block.sum$mb
  mb_ctrl = block.sum$mb_ctrl
  units.block = block.sum$units.block
  block.levels = block.sum$block.levels
  
  if(is.null(score.list.all)){
    score.list.all = score_all_blocks(nb, method.list.all)
  }
  
  result = 0
  for(i in 1 : B){
    score = score.list.all[[i]]
    
    
    zs = Z[units.block[[i]]]
    ys = Y[units.block[[i]]]
    
    stat.block = sum( score[rank(ys, ties.method = "first")[zs == 1] ] )
    result = result + weight[i] * 1 / mb[i] * stat.block
  }
  return(result)
}

#' function calculating mean and standard deviation
#' of each individual stratified rank sum statistic
#' 
#' output is a list of mean and standard deviations of each rank sum statistic

mu_sigma_single_weight_strat_rank_sum_stat <- function(Z, block,
                                                       method.list.all,
                                                       score.list.all = NULL,
                                                       weight,
                                                       block.sum = NULL){

  if(is.null(block.sum)){
    block.sum = summary_block(Z, block)
  }
  block = block.sum$block
  B = block.sum$B
  nb = block.sum$nb
  mb = block.sum$mb
  mb_ctrl = block.sum$mb_ctrl
  units.block = block.sum$units.block
  block.levels = block.sum$block.levels
  
  if(is.null(score.list.all)){
    score.list.all = score_all_blocks(nb, method.list.all)
  }

  
  mean_block = rep(NA, B)
  var_block = rep(NA, B)
  
  for(b in 1 : B){
    method.list = method.list.all[[b]]
    ns = nb[b]
    ms = mb[b]
    cs = mb_ctrl[b]
    
    score_b = score.list.all[[b]]
    mean_block[b] = mean(score_b)
    var_block[b] = cs/(ms*ns) * var(score_b)
  }
  
  single_mu = sum(mean_block * weight)
  
  single_sigma = sqrt(sum(weight^2 * var_block))
  return( list(mu = single_mu, sigma = single_sigma) )
}


#' generating null distribution of maximum among 
#' (1) weighted, (2) normalized multiple stratified rank sum statistic
#' 
#' output is a vector approximating distribution of null distribution of maximum among test statistics.
com_null_dist_block <- function(Z, block,
                                methods.list.all,
                                scores.list.all = NULL,
                                null.max = 10^4,
                                weight,
                                block.sum = NULL,
                                Z.perm = NULL){
  
  n = length(Z)
  H = length(methods.list.all)
  
  if(is.null(block.sum)){
    block.sum = summary_block(Z, block)
  }
  block = block.sum$block
  B = block.sum$B
  nb = block.sum$nb
  mb = block.sum$mb
  mb_ctrl = block.sum$mb_ctrl
  units.block = block.sum$units.block
  block.levels = block.sum$block.levels
  
  if(is.null(scores.list.all)){
    scores.list.all = list()
    for(h in 1:H){
      method.list.all = methods.list.all[[h]]
      scores.list.all[[h]] = score_all_blocks(nb, method.list.all)
    }
  }
  
  ## calculating mean and sigma of each stratified rank sum statistic
  mu_vec = rep(NA, H)
  sig_vec = rep(NA, H)
  for(h in 1:H){
    method.list.all = methods.list.all[[h]]
    score.list.all = scores.list.all[[h]]
    tmp = mu_sigma_single_weight_strat_rank_sum_stat(Z, block, method.list.all, score.list.all, weight, block.sum)
    
    mu_vec[h] = tmp$mu
    sig_vec[h] = tmp$sigma
  }
  
  ## generating random assignment vector corresponding to given block structure
  Z.perm = assign_block(block.sum, null.max)
  
  test.stat.all = matrix(NA, nrow = H, ncol = null.max)
  Y = c(1:n)
  
  for(iter2 in 1 : null.max){
    Z.tmp = Z.perm[, iter2]
    for(h in 1 : H){
      method.list.all = methods.list.all[[h]]
      tmp2 = single_weight_strat_rank_sum_stat(Z = Z.tmp, Y,
                                               block = block,
                                               method.list.all = method.list.all,
                                               score.list.all = scores.list.all[[h]],
                                               weight = weight,
                                               block.sum = block.sum)
      
      test.stat.all[h ,iter2] = (tmp2 - mu_vec[h]) / sig_vec[h]
    }
  }
  
  test.stat.max = apply(test.stat.all, 2, max)
  return(test.stat.max)
}


#' function for calculating population mean / standard deviation for a single
#' rank sum statistics
#' XL: I change weight to a vector

mu_sigma_list = function(Z, block, weight, methods.list.all, scores.list.all = NULL, block.sum = NULL){
  if(is.null(block.sum)){
    block.sum = summary_block(Z, block)
  }
  block = block.sum$block
  B = block.sum$B
  nb = block.sum$nb
  mb = block.sum$mb
  mb_ctrl = block.sum$mb_ctrl
  units.block = block.sum$units.block
  block.levels = block.sum$block.levels
  
  
  H = length(methods.list.all)
  mu = rep(0, H)
  sigma = rep(0, H)
  
  
  if(is.null(scores.list.all)){
    scores.list.all = list()
    for(h in 1:H){
      method.list.all = methods.list.all[[h]]
      scores.list.all[[h]] = score_all_blocks(nb, method.list.all)
    }
  }
  
  for(h in 1:H){
    method.list.all = methods.list.all[[h]]
    score.list.all = scores.list.all[[h]]
    temp = mu_sigma_single_weight_strat_rank_sum_stat(Z, block, method.list.all, score.list.all, weight, block.sum)
    mu[h] = temp$mu
    sigma[h] = temp$sigma
  }
  return(list(mu = mu, sigma = sigma))
}



#' calculating score values for given multiple stratified rank sum statistics
#' 
#' output is a list with H elements, 
#' each of the H elements is also a list of length B, 
#' which is a 2*(1+mb) matrix containing the minimum value of the block-specific 
#' test statistic given 0:mb numbers of units with effects greater than c

comb_matrix_block <- function(Z, Y, block, c, methods.list.all, scores.list.all = NULL, block.sum = NULL){
  if(is.null(block.sum)){
    block.sum = summary_block(Z, block)
  }
  block = block.sum$block
  B = block.sum$B
  nb = block.sum$nb
  mb = block.sum$mb
  mb_ctrl = block.sum$mb_ctrl
  units.block = block.sum$units.block
  block.levels = block.sum$block.levels
  

  H = length(methods.list.all)
  
  N = length(Y)
  
  if(is.null(scores.list.all)){
    scores.list.all = list()
    for(h in 1:H){
      method.list.all = methods.list.all[[h]]
      scores.list.all[[h]] = score_all_blocks(nb, method.list.all)
    }
  }
  
  total.list = list()
  for (l in 1 : H){
    method.list.all = methods.list.all[[l]]
    score.list.all = scores.list.all[[l]]
    Tlist = list()
    for(i in 1:B){
      Zb = Z[ units.block[[i]] ]
      Yb = Y[ units.block[[i]] ]
      Ti = matrix(nrow = 2, ncol = mb[i] + 1)
      Ti[1,] = 0:mb[i]
      
      for(ii in 0 : mb[i]){
        
        method.list = method.list.all[[i]]
        score = score.list.all[[i]]
        
        Ti[2, ii + 1] = min_stat(Zb, Yb, nb[i] - ii, c, method.list = method.list, score = score)
      }
      Tlist[[i]] = Ti
    }
    total.list[[l]] <- Tlist
  }
  return(total.list)
}


#' optimization function using integer programming by Gurobi
#' p is upper bound on the number of units with effects greater than c

Gurobi_sol_com <- function(Z, block, weight, coeflists, p, ms_list, exact = TRUE, block.sum = NULL){
  
  require(gurobi)
  require(Matrix)
  
  ## nn: number of all units + B
  ## B: number of strata
  ## H: number of stat to combine
  ## Q is organized by x's(n*H), eta's(H), theta
  ## weight is a H by B matrix, where each hth row is a weight vector for hth statistic.
  
  model = list()
  H = length(coeflists)
  B = length(coeflists[[1]])
  
  if(is.null(block.sum)){
    block.sum = summary_block(Z, block)
  }
  block = block.sum$block
  mb = block.sum$mb
  mb_ctrl = block.sum$mb_ctrl
  units.block = block.sum$units.block
  block.levels = block.sum$block.levels
  
  mmb = mb + 1 ## 0 ~ m_s + 1; for the constraint of sum of the number of k's in each strata
  
  
  ## parameter setting in the gurobi function
  
  nn = sum(mmb) 
  Q = rep(0, nn + H + 1)
  Q[nn + H + 1] = 1
  
  ## sum of x in each strata is 1
  Ai = rep(1:B, mmb)
  Aj = 1:nn
  x = rep(1, nn)
  
  # cost is p for all choices of stat and profit for h's stat is eta'h
  
  ## choosing x_{sj} for each stratum s, where x_{sj}=1 if and only if j units in this stratum have effects greater than c
  Ai = c(Ai, rep((B + 1), nn))
  Aj = c(Aj, (1:nn) )
  coeflist = coeflists[[1]]
  for(i in 1 : B){
    x = c(x, coeflist[[i]][1,])
  }
  
  ## listing t_s^h(j)s
  for (h in 1 : H){
    coeflist = coeflists[[h]]
    Ai = c(Ai, rep(B+1+h, nn))
    Aj = c(Aj, (1:nn) )
    for (b in 1 : B){
      x = c(x, weight[b] * (1 / mb[b]) * coeflist[[b]][2,])
    }
  }
  
  
  Ai = c(Ai, (B+1+1): (B+1+H))
  Aj = c(Aj, (nn+1) : (nn + H))
  x = c(x, rep(-1, H))
  
  mu = ms_list$mu
  sig_rev = 1/ms_list$sigma
  
  ## for final constraints
  
  Ai = c(Ai, rep((B + H + 1 + 1): (B + 2 * H + 1) , 2 ) )
  Aj = c(Aj, rep(nn + H + 1, H), (nn+1): (nn+H))
  x = c(x, rep(-1, H), sig_rev)
  
  A = Matrix::sparseMatrix(Ai, Aj, x = x)
  
  
  model$A = A
  model$obj = Q
  model$modelsense = "min"
  model$rhs = c(rep(1, B), p, rep(0, H), mu*sig_rev)
  model$sense = c(rep("=", B), "<=", rep("=", H), rep("<=", H))
  model$lb = c(rep(0, nn), rep(-Inf, H + 1))
  if(exact){
    model$vtype = c(rep("B", nn), rep("C", H + 1 ))
  }
  
  
  #  gurobi(model)
  params <- list(OutputFlag = 0)
  result = gurobi::gurobi(model,params)
  
  return(list(sol = result$x, obj = result$objval ))
}



#' Combined test statistic and combined p-value for randomization test
#' for quantiles of individual treatment effects
#' 
#' Obtain the combined stratified rank sum statistic
#' and combined p-value for testing the given
#' null hypothesis H0: \eqn{\tau_{(k)} \leq c}, or \eqn{H0: \tau_{(k)} \geq c}, or
#' \eqn{H0: \tau_{(k)} = c}, where \eqn{\tau_{(k)}} denotes individual
#' treatment effect at rank k.
#' 
#' @param Z An \eqn{n} dimensional treatment assignment vector.
#' @param Y An \eqn{n} dimensional observed outcome vector.
#' @param k An integer between 1 and n specifying which quantile of
#'   individual effect is of interest.
#' @param c A numerical object specifying the threshold for the null
#'   hypothesis.
#' @param block An \eqn{n} dimensional vector specifying block of each
#'   units. 
#' @param methods.list.all A list of lists of lists. Corresponds to the
#'.  method for each stratum for each different stratified rank sum 
#'.  statistic.
#' @param weight.name Weighting method to be implemented. If
#'   "asymp.opt", asymptotically optimal scheme under a class of local 
#'   alternatives is adjusted, if "dist.free", design-free scheme is 
#'   adjusted.
#' @param stat.null An vector whose empirical distribution
#'   approximates the randomization distribution of the combined 
#'   stratified rank sum statistic.
#' @param null.max A positive integer representing the number of
#'   permutations for approximating the randomization distribution of
#'   the rank sum statistic.
#' @param statistic logical; if TRUE (default), also prints the 
#'   combined stratified rank sum statistic.
#' @param opt.method Optimization method. If "ILP_gurobi", implement
#'   integer linear programming. If "LP_gurobi", implement relaxed
#'   linear programming. 
#'  
#' @return The p-value (and test statistic) for testing the specified
#'   null hypothesis of interest.
#' @export  
pval_comb_block <- function(Z, Y, k, c,
                            block,
                            methods.list.all,
                            weight.name = "asymp.opt",
                            stat.null = NULL,
                            null.max = 10^5,
                            statistic = TRUE,
                            opt.method = "ILP_gurobi"){
  
  if(opt.method == "ILP_gurobi"){
    exact = TRUE
  }
  
  if(opt.method == "LP_gurobi"){
    exact = FALSE
  }
  
  H = length(methods.list.all)
  
  N = length(Y)
  p = N - k
  
  block.sum = summary_block(Z, block)
  block = block.sum$block
  
  weight = weight_scheme(block.sum, weight.name)
  scores.list.all = list()
  for(h in 1:H){
    method.list.all = methods.list.all[[h]]
    scores.list.all[[h]] = score_all_blocks(block.sum$nb, method.list.all)
  }
  
  ms_list = mu_sigma_list(Z, block, weight, methods.list.all, scores.list.all, block.sum)
  
  if(is.null(stat.null)){
    stat.null = com_null_dist_block(Z, block, methods.list.all, scores.list.all, null.max, weight, block.sum, Z.perm = NULL)
  }
  
  coeflists = comb_matrix_block(Z, Y, block, c, methods.list.all, scores.list.all, block.sum)
  
  
  stat.min = Gurobi_sol_com(Z, block, weight, coeflists, p, ms_list, exact, block.sum)$obj
  
  pval = mean(stat.null >= stat.min)
  
  if(statistic == TRUE) {
    result = c(pval, stat.min)
    names(result) = c("p.value", "test.stat")
    return(result)
  } else return(pval)
}



#' Helper function for Simultaneous Inference for multiple quantiles on CRE
#' 
#' Output is a lower limit of prediction intervals for prespecified quantiles.

com_block_conf_quant_larger_trt <- function(Z, Y,
                                            block,
                                            k.vec = NULL,
                                            methods.list.all,
                                            weight.name = "asymp.opt",
                                            opt.method = "ILP_gurobi",
                                            stat.null = NULL, null.max = 10^4,
                                            tol = 0.01,
                                            alpha = 0.1){
  
  if(opt.method == "ILP_gurobi"){
    exact = TRUE
  }
  
  if(opt.method == "LP_gurobi"){
    exact = FALSE
  }
  
  n = length(Z)
  m = sum(Z)
  
  block.sum = summary_block(Z, block)
  block = block.sum$block
  B = block.sum$B
  
  H = length(methods.list.all)
  
  weight = weight_scheme(block.sum, weight.name)
  scores.list.all = list()
  for(h in 1:H){
    method.list.all = methods.list.all[[h]]
    scores.list.all[[h]] = score_all_blocks(block.sum$nb, method.list.all)
  }
  
  ms_list = mu_sigma_list(Z, block, weight, methods.list.all, scores.list.all, block.sum)

  if(is.null(stat.null)){
    stat.null = com_null_dist_block(Z, block, methods.list.all, scores.list.all, null.max, weight, block.sum, Z.perm = NULL)
  }
  
  thres = sort(stat.null, decreasing = TRUE)[floor(null.max * alpha) + 1]
  
  
  Y1.max = max(Y[Z==1])
  Y1.min = min(Y[Z==1])
  Y0.max = max(Y[Z==0])
  Y0.min = min(Y[Z==0])
  
  c.max = Y1.max - Y0.min + tol
  c.min = Y1.min - Y0.max - tol
  
  if( is.null(k.vec) ){
    
    quantiles = rep(NA, n)
    
    for(k in n : (n-m)){
      
      if(k < n){
        c.max = quantiles[k+1]
      }
      p = n - k
      
      f <- function(c){
        coeflists = comb_matrix_block(Z, Y, block, c, methods.list.all, scores.list.all, block.sum)
        stat.min = Gurobi_sol_com(Z, block, weight, coeflists, p, ms_list, exact, block.sum)$obj
        return(stat.min - thres)
      }
      
      tmp1 = f(c.min)
      tmp2 = f(c.max)
      if( tmp1 <= 0){
        c.sol = -Inf
      }
      if( tmp2 > 0){
        c.sol = c.max
      }
      if( tmp1 > 0 & tmp2 <= 0 ){  # maybe possible to improve this by using connecting previous results logically
        c.sol = uniroot(f, interval = c(c.min, c.max), extendInt = "downX", tol = tol)$root
        c.sol = round(c.sol, digits = -log10(tol))
        
        if( f(c.sol) <= 0 ){ #### need to update this, following Professor's changes on CRE
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
      
      if( c.sol == -Inf & k > (n-m)) {
        quantiles[ (n-m) : (k-1)] = -Inf
        break
      }
    }
    
    if( n-m > 1){
      quantiles[1 : (n-m-1)] = quantiles[n-m]
    }
    
    quantiles[quantiles > (Y1.max - Y0.min) + tol/2] = Inf
    
  }else{
    
    k.vec.sort = sort(k.vec, decreasing = FALSE)
    j.max = length(k.vec.sort)
    j.min = max( sum(k.vec <= (n-m) ), 1)
    
    quantiles = rep(NA, j.max)
    
    for(j in j.max:j.min){
      
      k = k.vec.sort[j]
      p = n - k
      if(j < j.max){
        c.max = quantiles[j+1]
      }
      
      f <- function(c){
        coeflists = comb_matrix_block(Z, Y, block, c, methods.list.all, scores.list.all, block.sum)
        stat.min = Gurobi_sol_com(Z, block, weight, coeflists, p, ms_list, exact, block.sum)$obj
        return(stat.min - thres)
      }
      
      tmp1 = f(c.min)
      tmp2 = f(c.max)
      if( tmp1 <= 0){
        c.sol = -Inf
      }
      if( tmp2 > 0){
        c.sol = c.max
      }
      if( f(c.min) > 0 & f(c.max) <= 0 ){
        c.sol = uniroot(f, interval = c(c.min, c.max), extendInt = "downX",
                        tol = tol)$root
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
      
      if (c.sol == -Inf & j > j.min){
        quantiles[j.min : (j-1)] = -Inf
        break
      }
    }
    
    if(j.min > 1){
      quantiles[1 : (j.min - 1)] = quantiles[j.min]
    }
    
  }
  return(quantiles)
}


#' Simultaneous bound for confidence interval using combined rank sum statistic 
#' on stratified randomized experiment
#' 
#' null hypothesis H0: \eqn{\tau_{(k)} \leq c}, or \eqn{H0: \tau_{(k)} \geq c}, or
#' \eqn{H0: \tau_{(k)} = c}, where \eqn{\tau_{(k)}} denotes individual
#' treatment effect at rank k.
#' 
#' @param Z An \eqn{n} dimensional treatment assignment vector.
#' @param Y An \eqn{n} dimensional observed outcome vector.
#' @param k An integer between 1 and n specifying which quantile of
#'   individual effect is of interest.
#' @param c A numerical object specifying the threshold for the null
#'   hypothesis.
#' @param block An \eqn{n} dimensional vector specifying block of each
#'   units. 
#' @param set set of quantiles of interests. If "treat", generate prediction intervals
#'   for effect quantiles among treated units. If "control", generate prediction intervals
#'   for effect quantiles among control units. If "all", generate confidence intervals for
#'   all effect quantiles.
#' @param methods.list.all A list of lists of lists. Corresponds to the
#'.  method for each stratum for each different stratified rank sum 
#'.  statistic.
#' @param weight.name Weighting method to be implemented. If
#'   "asymp.opt", asymptotically optimal scheme under a class of local 
#'   alternatives is adjusted, if "dist.free", design-free scheme is 
#'   adjusted.
#' @param opt.method Optimization method. If "ILP_gurobi", implement
#'   integer linear programming. If "LP_gurobi", implement relaxed
#'   linear programming. 
#' @param stat.null An vector whose empirical distribution
#'   approximates the randomization distribution of the combined 
#'   stratified rank sum statistic.
#' @param null.max A positive integer representing the number of
#'   permutations for approximating the randomization distribution of
#'   the rank sum statistic.
#' @param tol A numerical object specifying the precision of the
#'   obtained confidence intervals. For example, if tol = 10^(-3),
#'   then the confidence limits are precise up to 3 digits.
#' @param alpha A numerical object, where 1-alpha indicates the
#'   confidence level.
#' @return A vector specifying lower limits of prediction (confidence) intervals for
#' quantiles k = 1 ~ m (or n - m, or n).


com_block_conf_quant_larger = function(Z, Y,
                                       block,
                                       set = "all",
                                       methods.list.all = NULL,
                                       weight.name = "asymp.opt",
                                       opt.method = "ILP_gurobi",
                                       stat.null = NULL, 
                                       null.max = 10^4, 
                                       tol = 0.01, 
                                       alpha = 0.1){
  
  n = length(Z)
  
  if(set == "treat"){
    ci.treat = com_block_conf_quant_larger_trt(Z, Y, block,
                                               k.vec = NULL,
                                               methods.list.all,
                                               weight.name,
                                               opt.method,
                                               stat.null, null.max, tol, alpha)
    ci.treat = ci.treat[(n - sum(Z) + 1) : n]
    return(ci.treat)
  }
  
  if(set == "control"){
    Z = 1 - Z
    Y = -Y
    ci.control = com_block_conf_quant_larger_trt(Z, Y, block,
                                                 k.vec = NULL,
                                                 methods.list.all,
                                                 weight.name,
                                                 opt.method,
                                                 stat.null, null.max, tol,
                                                 alpha)
    ci.control = ci.control[(n - sum(Z) + 1) : n]
    return(ci.control)
  }
  
  if(set == "all"){
    ci.treat = com_block_conf_quant_larger_trt(Z, Y, block,
                                               k.vec = NULL,
                                               methods.list.all,
                                               weight.name,
                                               opt.method,
                                               stat.null, null.max, tol,
                                               alpha = alpha / 2)
    ci.treat = ci.treat[(n - sum(Z) + 1) : n]
    
    Z = 1 - Z
    Y = -Y
    
    ci.control = com_block_conf_quant_larger_trt(Z, Y, block,
                                                 k.vec = NULL,
                                                 methods.list.all,
                                                 weight.name,
                                                 opt.method,
                                                 stat.null, null.max, tol,
                                                 alpha = alpha / 2)
    ci.control = ci.control[(n - sum(Z) + 1) : n]
    
    ci.all = sort( c(ci.treat, ci.control))
    return(ci.all)
  }
}


