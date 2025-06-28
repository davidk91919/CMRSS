###################################
## helper functions
###################################

#sorting index among treated units, according to ranks of the outcome
#output: indices of treated units with increasing outcomes
sort_treat <- function(Y, Z){
  r = rank(Y, ties.method = "first")
  ind.sort = sort.int(r, index.return = TRUE)$ix
  ind.sort.treat = ind.sort[Z[ind.sort] == 1]
  return(ind.sort.treat)
}

#Generate matrix of complete randomization assignments
#output: a matrix with nperm columns, each column gives one assignment
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


# calculating score for given n
# included polynomial rank score with Puri(1965)'s normalization (denoted by std)
# output: rank scores of length n
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

# generating randomization null distribution of a single rank sum statistic
# output: null dist of a single rank score statistics under the CRE
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

# XL: I add the following
# null distribution of multiple rank sum statistics
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


# calculating individual test statistic corresponds to Caughey et al.(2023)
# output: minimum value of the rank sum statistic assuming n-k treated units can have infinite effects
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

# calculating p-value for each individual methods
#output: the p-value assuming at most n-k treated units having effects > c, using a single rank sum statistic
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

###################################
## CRE with multiple rank statistics
###################################

# null distribution for minimum p-values
# output: the distribution of the minimum p-value, particularly Monte Carlo samples from \overline{F} in Theorem 1 of the paper
comb_null_dist_cre = function(n, m, methods.list, Z.perm = NULL, nperm = 10^4, stat.null.mult = NULL){

  H = length(methods.list)

  if(is.null(stat.null.mult)){
    if(is.null(Z.perm)){
      Z.perm = assign_CRE(n, m, nperm)
    }
    stat.null.mult = null_dist_multiple(n, m, methods.list, nperm, Z.perm)
  }

  # nperm = ncol(Z.perm)
  # stat.dist = matrix(NA, nrow = H, ncol = nperm)
  tail.prob = matrix(NA, nrow = H, ncol = nperm)
  #
  # for(i in 1 : H){
  #   method.list = methods.list[[i]]
  #   stat.dist[i,] = null_dist(n, m, method.list, Z.perm = Z.perm)
  # }

  for(j in 1 : H){
#    prob.dist[j,] = rank(stat.dist[j,]) / nperm
#    prob.dist[j,] = 1 -  (rank(stat.dist[j,]) - 1) / nperm ############
    # for(l in 1 : nperm){
    #   tail.prob[j, l] = mean( stat.dist[j,] >= stat.dist[j,l] ) # double check whether they are the same
    # }
    tail.prob[j, ] = 1 - (rank(stat.null.mult[j,], ties.method = "min") - 1) / nperm
  }

  result = apply(tail.prob, 2, min)
  return(result)
}


### calculate the minimum p value from multiple rank sum statistics
min_p_multiple_rank_sum <- function(Z, Y, k, c, methods.list, Z.perm = NULL, nperm, stat.null.mult = NULL){
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

# main function1 : function for calculating valid minimum p-value, corresponds to Corollary 1 in the draft
# output: a valid p-value, based on multiple random sum statistics, for testing the null that at most n-k treated units having effects > c
# I change to combined p value, to distinguish it from min p value
comb_p_val_cre = function(Z, Y, k, c, methods.list, stat.null = NULL, Z.perm = NULL, nperm,
                             stat.null.mult = NULL,
                             null.dist.comb = NULL
                             ){ # XL: I add Z.perm in the input and reorder the input

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

  # if(is.null(stat.null)){
  #   null.dist = comb_null_dist_cre(n = n, m = m,
  #                                  methods.list = methods.list, Z.perm = Z.perm, nperm = nperm)
  # }else{
  #   null.dist = stat.null
  # }

  # XL: stat.null is fine, no need for null.dist. I also change the "result" later
  # stat.null gives Monte Carlo samples from the null distribution of the minimum p-value
  # if(is.null(stat.null)){
  #   stat.null = comb_null_dist_cre(n = n, m = m,
  #                                  methods.list = methods.list, Z.perm = Z.perm, nperm = nperm)
  # }
  # # p-value for each single rank sum statistic
  # for (i in 1 : H){
  #   method.list = methods.list[[i]]
  #   pval.vec[i] = pval_cre(Z = Z, Y = Y, k = k, c = c, Z.perm = Z.perm, method.list = method.list)
  # }
  #
  # min.pval = min(pval.vec)

  min.pval = min_p_multiple_rank_sum(Z, Y, k, c, methods.list, Z.perm, nperm, stat.null.mult)

  # result = mean(null.dist <= min.pval)
  result = mean(null.dist.comb <= min.pval) # XL: I changed
  return(result)
}

# function for Simultaneous Inference for multiple quantiles on CRE
com_conf_quant_larger_trt <- function( Z, Y, methods.list = NULL,
                                       # stat.null = NULL,
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
  ## DK: Is this just because to get matched sign in the future part?
  thres.minus.min.pval = sort(-1*null.dist.comb, decreasing = TRUE)[ floor(nperm * alpha) + 1 ]

  #XL??? is this the right threshold? Should it be
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

      ## the original way
    #   f <- function(c){
    #     pval = comb_p_val_cre(Z, Y, k, c,  methods.list = methods.list, stat.null = stat.null, Z.perm = Z.perm, nperm = nperm) # not a min_stat, but a min pval.
    #     #XL: it is better to call it a combined p-value, distinguishing it from min pval
    #     return(thres - pval)
    #   }
    #   # check whether f(-Inf) = f(c.min) > 0 #
    #   if( f(c.min) <= 0 ){  ##
    #     c.sol = -Inf
    #   }
    #   if( f(c.max) > 0){  ##
    #     c.sol = c.max
    #   }
    #   if( f(c.min) > 0 & f(c.max) <= 0 ){   # need to change this in reverse (c.min, c.max)
    #     c.sol = uniroot(f, interval = c(c.min, c.max), extendInt = "downX", tol = tol)$root
    #     # find the min c st p-value > alpha <==> f <= 0 #
    #     c.sol = round(c.sol, digits = -log10(tol))
    #     # if( f(c.sol) <= 0 ){
    #     #   while( f(c.sol) <= 0 ){
    #     #     c.sol = c.sol - tol
    #     #   }
    #     #   c.sol = c.sol + tol
    #     # }else{
    #     #   while(f(c.sol) > 0 ){
    #     #     c.sol = c.sol + tol
    #     #   }
    #     # }
    #     ###XL??? I change to the following
    #     if( f(c.sol) < 0 ){
    #       while( f(c.sol) < 0 ){
    #         c.sol = c.sol - tol
    #       }
    #       c.sol = c.sol + tol
    #     }else{
    #       while(f(c.sol) >= 0 ){
    #         c.sol = c.sol + tol
    #       }
    #     }
    #   }
      c.limit[k] = c.sol
    }

    if( n-m > 1){
      c.limit[1:(n-m-1)] = c.limit[n-m]
    }

    #c.limit[c.limit > (Y1.max - Y0.min) + tol/2] = Inf
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

      # the original way
      # f = function(c){
      #   pval = min_p_val_cre(Z, Y, k, c, stat.null, methods.list, nperm)
      #   return(thres - pval)
      # }
      # if( f(c.min) <= 0){
      #   c.sol = -Inf
      # }
      # if( f(c.max) > 0){
      #   c.sol = c.max
      # }
      # if(f(c.min) > 0 & f(c.max) <= 0){
      #   c.sol = uniroot(f, interval = c(c.min, c.max), extendInt = "downX", tol = tol)$root
      #   c.sol = round(c.sol, digits = -log10(tol))
      #   # if( f(c.sol) <= 0){
      #   #   while( f(c.sol) <= 0){
      #   #     c.sol = c.sol - tol
      #   #   }
      #   #   csol = c.sol + tol
      #   # }else{
      #   #   while(f(c.sol) > 0 ){
      #   #     c.sol = c.sol + tol
      #   #   }
      #   # }
      #   # XL???: I change to the following
      #   if( f(c.sol) < 0){
      #     while( f(c.sol) < 0){
      #       c.sol = c.sol - tol
      #     }
      #     csol = c.sol + tol
      #   }else{
      #     while(f(c.sol) >= 0 ){
      #       c.sol = c.sol + tol
      #     }
      #   }
      #
      # }
      c.limit[j] = c.sol
    }
    if(j.min > 1){
      c.limit[1:(j.min-1)] = c.limit[j.min]
    }
    #    c.limit[c.limit > (Y1.max - Y0.min) + 10^{-log10(tol)} / 2] = Inf
  }

  return( c.limit )
}

# main function2 : Simultaneous lower bounds using multiple rank sum statistics in CRE
# XL: I delete some input of the function
com_conf_quant_larger_cre <- function( Z, Y, methods.list,
                                      # stat.null = NULL, # by now, this is just for the treated ones(will add option for scores for control units)
                                      nperm = 10^4,
                                      set = "treat",
                                      Z.perm = NULL,
                                      # k.vec = NULL, ## DK: don't we need an option for the specific quantiles?
                                      alpha = 0.05,
                                      tol = 10^(-3)
                                      # ind.sort.treat = NULL
                                      ){

  H = length(methods.list)
  n = length(Z)
  if(is.null(Z.perm)){
    Z.perm = assign_CRE(n, m=sum(Z), nperm)
  }
  # set = "prev": without Chen and Li's adjustment
  # if(set == "prev"){
  #   ci.prev = com_conf_quant_larger_trt(Z, Y, methods.list, stat.null,
  #                                       nperm, k.vec,
  #                                       Z.perm, alpha, tol, ind.sort.treat
  #   )
  #   return(ci.prev)
  # }

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


###################################
## SERE; from individual p-values to simultaneous lower bound for multiple quantiles
###################################


# XL: I add some basic functions
### summary of the block
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

### give the weight under given scheme
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

### score for each of the blocks under a given method.list.all
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

# Workflow1 : Generating null distribution in SERE

## helper1: function to calculate t in (14) of the paper for a single stratified rank sum (without mu, sigma involved)

# XL: I define weight directly as a vector below

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

  # if(!is.factor(block)){
  #   block = as.factor(block)
  # }
  # block.levels = levels(block)
  # B = length(block.levels)
  #
  # # nb = vector(length = B)
  # #XL??-a minor change
  # nb = rep(NA, B)
  #
  # mb = rep(NA, B)
  #
  # mb_ctrl = rep(NA, B)
  # for(i in 1 : B){
  #   nb[i] = length(Z[block == block.levels[i]])
  #   mb[i] = sum(Z[block == block.levels[i]])
  # }
  # mb_ctrl = nb - mb

  # if(weight == "asymp.opt"){
  #   weight = mb
  # } else{
  #   if(weight == "dis.free"){
  #     weight = (nb + 1) / mb_ctrl
  #   }
  # }

  result = 0
  for(i in 1 : B){
    # method.list = method.list.all[[i]]
    # score = rank_score(nb[i], method.list)
    score = score.list.all[[i]]

    # zs = Z[block == block.levels[[i]]]
    # ys = Y[block == block.levels[[i]]]

    zs = Z[units.block[[i]]]
    ys = Y[units.block[[i]]]

    stat.block = sum( score[rank(ys, ties.method = "first")[zs == 1] ] )
    result = result + weight[i] * 1 / mb[i] * stat.block
  }
  return(result)
}

## helper2: function calculating mu and sigma of each individual stratified rank sum statistic
# XL: change weight directly to a vector
mu_sigma_single_weight_strat_rank_sum_stat <- function(Z, block, method.list.all, score.list.all = NULL, weight, block.sum = NULL){
  #################################
  ######## started from here ######
  #################################

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

  # if(!is.factor(block)){
  #   block = as.factor(block)
  # }
  # block.levels = levels(block)
  # B = length(block.levels)
  #
  # nb = vector(length = B)
  # mb = rep(NA, B)
  #
  # mb_ctrl = rep(NA, B)
  # for(i in 1 : B){
  #   nb[i] = length(Z[block == block.levels[i]])
  #   mb[i] = sum(Z[block == block.levels[i]])
  # }
  # mb_ctrl = nb - mb

  # if(weight == "asymp.opt"){
  #   weight = mb
  # } else{
  #   if(weight == "dis.free"){
  #     weight = (nb + 1) / mb_ctrl
  #   }
  # }

  mean_block = rep(NA, B)
  var_block = rep(NA, B)

  for(b in 1 : B){
    method.list = method.list.all[[b]]
    ns = nb[b]
    ms = mb[b]
    cs = mb_ctrl[b]

    # score_b = rank_score(ns, method.list)
    score_b = score.list.all[[b]]

    # mean_block[b] = 1 / ns * sum(score_b)
    # var_block[b] = cs/(ms*ns*(ns-1)) * sum((score_b - mean_block[b])^2)

    mean_block[b] = mean(score_b)
    var_block[b] = cs/(ms*ns) * var(score_b)
  }

  single_mu = sum(mean_block * weight)

  ##XL??? could put sqrt() here instead of in return
  single_sigma = sqrt(sum(weight^2 * var_block))

  return( list(mu = single_mu, sigma = single_sigma) )
}

## generating null distribution of maximum among (1) weighted, (2) normalized multiple stratified rank sum statistic
# XL: I change weight to a vector
com_null_dist_block <- function(Z, block, methods.list.all, scores.list.all = NULL, null.max = 10^4, weight, block.sum = NULL, Z.perm = NULL){

  # if(!is.factor(block)){
  #   block = as.factor(block)
  # }
  # block.levels = levels(block)
  # B = length(block.levels)
  #
  # n = length(Z)
  # nb = vector(length = B)
  # mb = rep(NA, B)
  #
  # for(i in 1 : B){
  #   nb[i] = length(Z[block == block.levels[i]])
  #   mb[i] = sum(Z[block == block.levels[i]])
  # }

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
    # tmp = mu_sigma_single_weight_strat_rank_sum_stat(Z = Z, block = block, method.list.all = method.list.all, weight = weight)
    tmp = mu_sigma_single_weight_strat_rank_sum_stat(Z, block, method.list.all, score.list.all, weight, block.sum)

    mu_vec[h] = tmp$mu
    sig_vec[h] = tmp$sigma
  }

  ## generating random assignment vector corresponding to given block structure
  Z.perm = assign_block(block.sum, null.max)

  # Z.perm = matrix(NA, nrow = n, ncol = null.max)
  # for(iter in 1 : null.max){
  #   tmp = rep(0, n)
  #   for(i in 1 : B){
  #     tmp[which(block == block.levels[i])] = sample( c(rep(1, mb[i]), rep(0, nb[i] - mb[i])) )
  #   }
  #   Z.perm[,iter] = tmp
  # }

  ##
  test.stat.all = matrix(NA, nrow = H, ncol = null.max)
  Y = c(1:n)
  # rnorm(n) # getting any fixed outcome vector. just for giving input for the single test statistic; and will no matter since our test statistic is still distribution free.

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

  # check whether the null dist for each single stratified statistic has mean 0 and variance 1
  # apply(test.stat.all, 1, mean)
  # apply(test.stat.all, 1, sd)

  test.stat.max = apply(test.stat.all, 2, max)
  return(test.stat.max)
}

########################################################

# Workflow2 : Calculating minimum test statistic

## function for calculating population mean / variance for each statistics
#XL??? Why not directly use mu_sigma_single_weight_strat_rank_sum_stat???
# XL: I change weight to a vector
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

  # if(!is.factor(block)){
  #   block = as.factor(block)
  # }

  H = length(methods.list.all)
  mu = rep(0, H)
  sigma = rep(0, H)

  #XL??? define B using the length of block.levels, since B here represents the number of blocks
  #XL??? also, we can avoid s by using only B
  #  B = s = length(methods.list.all[[1]][[1]])
  # old one
  # B = s = length(methods.list.all[[1]])

  # block.levels = levels(block)
  # XL--changed
  # B = length(block.levels)
  # s = B

  # nb = rep(NA, B)
  # for (i in 1 :B){
  #   nb[i] = sum(block == block.levels[i])
  # }

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

  # > mu
  # [1] 3.023091 1.822342
  # > sigma
  # [1] 0.5497698 0.5482509

  # > mu
  # [1] 1.2228101 0.8046002
  # > sigma
  # [1] 0.2357893 0.2435407

# mean = rep(NA, H)
# if(weight == "asymp.opt"){
#   for(i in 1 : H){
#     method.list.all = methods.list.all[[i]]
#     tmp1 <- 0
#     tmp2 <- 0
#     for(j in 1 : B){
#       method.list = method.list.all[[j]]
#       ns = nb[j]
#       ms = sum(Z[block == block.levels[j]])
#
#       score = rank_score(ns, method.list)
#       tmp1 = tmp1 + ms * 1/ns * sum(score)
#       tmp2 = tmp2 + ms^2 * (1/ms - 1/ns) * (1 / (ns - 1)) *
#         sum( (score - 1 / ns * sum(score))^2 )
#       #       + (sum(score^2) - 1 / ns * sum(score)^2)      #same
#     }
#     mean[i] = tmp1
#     sigma[i] = sqrt(tmp2)
#     }
# }
#
# if(weight == "dis.free"){
#   for(i in 1 : H){
#     method.list.all = methods.list.all[[i]]
#     tmp1 <- 0
#     tmp2 <- 0
#     for(j in 1 : s){
#       method.list = method.list.all[[j]]
#       ns = nb[j]
#       ### XL??? this is not right, ms is the number of control, and you use ms as number of treated later
#       ###old
#       # ms = sum(1 - Z[block == block.levels[j]])
#       # w = (ns + 1) / ms
#       ###XL updated
#       ms = sum(Z[block == block.levels[j]])
#       ms_ctrl = sum(1 - Z[block == block.levels[j]])
#       w = (ns + 1) / ms_ctrl
#
#       score = rank_score(ns, method.list)
#       tmp1 = tmp1 + w * 1/ns * sum(score)
#       tmp2 = tmp2 + w^2 * (1/ms - 1/ns) * (1 / (ns - 1)) *
#         sum( (score - 1 / ns * sum(score))^2 )
#     }
#     mean[i] = tmp1
#     sigma[i] = sqrt(tmp2)
#   }
# }
  return(list(mu = mu, sigma = sigma))
}

## generating a list of lists, where each list corresponds to test statistics for each method(on each strata)
#XL--output: a list with H elements, each of the H elements is also a list of length B, which is a 2*(1+mb) matrix containing the minimum value of the block-specific test statistic given 0:mb numbers of units with effects greater than c

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


  # if(!is.factor(block)){
  #   block = as.factor(block)
  # }

  H = length(methods.list.all)

  N = length(Y)
  # B = length(levels(block))
  # block.levels = levels(block)

  ## calculate numbers of observations in each block
  # nb = rep(NA, B)
  # for(i in 1:B){
  #   nb[i] = sum(block == block.levels[i])
  # }

  ## calculate number of treated units in each block
  # mb = rep(NA, B)
  # for(i in 1 : B){
  #   mb[i] = sum(Z[block == block.levels[i]])
  # }

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
      #Ti[1,] = nb[i] - (0 : mb[i])
      for(ii in 0 : mb[i]){
        # if(length(method.list.all) == 1){
        #   method.list = method.list.all[[1]]
        # }else{

        method.list = method.list.all[[i]]
        score = score.list.all[[i]]

        # }
        Ti[2, ii + 1] = min_stat(Zb, Yb, nb[i] - ii, c, method.list = method.list, score = score)
      }
      Tlist[[i]] = Ti
    }
    total.list[[l]] <- Tlist
  }
  return(total.list)
}

#################################################################

# main function1 : optimization function using integer programming by Gurobi
# p is upper bound on the number of units with effects greater than c
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

  ## calculate numbers of treated units(m_s) / controlled units(mb_ctrl) in each block
  # mb = rep(NA, B)
  # for(i in 1:B){
  #   mb[i] = sum(Z[block == block.levels[i]])
  # }

  mmb = mb + 1 ## 0 ~ m_s + 1; for the constraint of sum of the number of k's in each strata

  # mb_ctrl = rep(NA, B)
  # for(i in 1:B){
  #   mb_ctrl[i] = sum(1 - (Z[block == block.levels[i]]))
  # }

  # weight matrix generation
  ## weight corresponds to scheme 1 ** default option
  # if(weight == "asymp.opt"){
  #   weight = matrix(nrow = H, ncol = B)
  #   for(i in 1 : H){
  #     weight[i,] = mb
  #   }
  # }
  # else{
  #   ## weight corresponds to scheme 2
  #   if(weight == "dis.free"){
  #     weight = matrix(nrow = H, ncol = B)
  #     for(i in 1 : H){
  #       weight[i,] = (mb + mb_ctrl + 1) / mb_ctrl
  #     }
  #   }
  # }

  ## parameter setting in the gurobi function

  nn = sum(mmb) ## considering 0 ~ m_s + 1; 0126 added
  #  indx = c(0, cumsum(rep(nb, length(coeflists))))
  Q = rep(0, nn + H + 1)
  Q[nn + H + 1] = 1
  # block.levels = levels(block)

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
      # no need for weight to be matrix
      # x = c(x, weight[h,b] * (1 / mb[b]) * coeflist[[b]][2,])
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
  # nnn+1 ~ nn+Hth column: etas / nn+H+1th column: t_*
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


## printing test statistic and p-value
# XL: change weight to weight.name
# I also change the function name
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

  N = length(Y)
  p = N - k

  block.sum = summary_block(Z, block)
  block = block.sum$block

  # if(!is.factor(block)){
  #   block = as.factor(block)
  # }

  weight = weight_scheme(block.sum, weight.name)
  scores.list.all = list()
  for(h in 1:H){
    method.list.all = methods.list.all[[h]]
    scores.list.all[[h]] = score_all_blocks(block.sum$nb, method.list.all)
  }

  ms_list = mu_sigma_list(Z, block, weight, methods.list.all, scores.list.all, block.sum)

  # ms_list = mu_sigma_list(Z = Z, block = block, weight = weight,
  #                         methods.list.all = methods.list.all)

  if(is.null(stat.null)){

    stat.null = com_null_dist_block(Z, block, methods.list.all, scores.list.all, null.max, weight, block.sum, Z.perm = NULL)

    # stat.null = com_null_dist_block(Z = Z, block = block,
    #                                 methods.list.all = methods.list.all,
    #                                 null.max = null.max
    #                                 #Z.perm.all = NULL,
    #                                 #ms.list = ms_list
    #                                 )
  }

  coeflists = comb_matrix_block(Z, Y, block, c, methods.list.all, scores.list.all, block.sum)

  # coeflists = comb_matrix_block(Z = Z, Y = Y,
  #                               block = block, c = c,
  #                               methods.list.all = methods.list.all)


  stat.min = Gurobi_sol_com(Z, block, weight, coeflists, p, ms_list, exact, block.sum)$obj

  # stat.min = Gurobi_sol_com(Z = Z, block = block,
                            # weight = weight,
                            # coeflists = coeflists,
                            # p = p,
                            # ms_list = ms_list,
                            # exact = exact)$obj

  pval = mean(stat.null >= stat.min)

  if(statistic == TRUE) {
    result = c(pval, stat.min)
    names(result) = c("p.value", "test.stat")
    return(result)
  } else return(pval)
}

##############################################################################


# Helper function for Simultaneous Inference for multiple quantiles on CRE

com_block_conf_quant_larger_trt <- function(Z, Y,
                                           block,
                                           #quantiles = NULL,
                                           k.vec = NULL,
                                           methods.list.all,
                                           weight.name = "asymp.opt",
                                           opt.method = "ILP_gurobi",
                                           stat.null = NULL, null.max = 10^4,
                                           tol = 0.01,
                                           alpha = 0.1){

  # if(!is.factor(block)){
  #   block = as.factor(block)
  # }
  # levels(block) = 1:length(levels(block))

  if(opt.method == "ILP_gurobi"){
    exact = TRUE
  }

  if(opt.method == "LP_gurobi"){
    exact = FALSE
  }

  n = length(Z)
  m = sum(Z)
  # B = length(levels(block))

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
  # ms.list = mu_sigma_list(Z, block, weight, methods.list.all)



  # if( all(is.nan(ms.list$mean)) ) { # XL: why we need this?
  #   stop("Check choice of test parameter is larger than the block size")
  # }

  if(is.null(stat.null)){
    stat.null = com_null_dist_block(Z, block, methods.list.all, scores.list.all, null.max, weight, block.sum, Z.perm = NULL)
  }

  thres = sort(stat.null, decreasing = TRUE)[floor(null.max * alpha) + 1]

  # Y1.max = max(Y[Z == 1])
  # Y0.min = min(Y[Z == 0])
  #
  # c.max = max(Y[Z == 1]) - min(Y[Z == 0]) + tol
  # c.min = min(Y[Z == 1]) - max(Y[Z == 0]) - tol

  Y1.max = max(Y[Z==1])
  Y1.min = min(Y[Z==1])
  Y0.max = max(Y[Z==0])
  Y0.min = min(Y[Z==0])

  c.max = Y1.max - Y0.min + tol
  c.min = Y1.min - Y0.max - tol

 #if( is.null(quantiles) ){
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
#      print(k) # for debugging. will be deleted in the final version
      # I add k > (n-m)
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
        # coeflists = comb_matrix_block(Z = Z, Y= Y,
        #                               block = block, c = c,
        #                               methods.list.all = methods.list.all)
        # stat.min = Gurobi_sol_com(Z = Z, block = block,
        #                           weight = weight,
        #                           coeflists = coeflists,
        #                           p = p,
        #                           ms_list = ms.list,
        #                           exact = exact)$obj
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
#      print(k) # for debugging
      # XL: I add j > j.min
      if (c.sol == -Inf & j > j.min){
        # XL: I edited, the original one is an error and can cause problem ## DK: maybe because different letter used in the loop?
        quantiles[j.min : (j-1)] = -Inf
        # quantiles[(n-m) : (k-1)] = -Inf
        break
      }
    }

    if(j.min > 1){
      quantiles[1 : (j.min - 1)] = quantiles[j.min]
    }

  }
  return(quantiles)
}

# main function 2: Simultaneous bound for confidence interval using CMRSS on SERE

com_block_conf_quant_larger = function(Z, Y,
                                       block,
                                       # k.vec = NULL,
                                       set = "all",
                                       methods.list.all = NULL,
                                       weight.name = "asymp.opt",
                                       opt.method = "ILP_gurobi",
                                       stat.null = NULL, null.max = 10^4, tol = 0.01, alpha = 0.1){

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

    # ci.trt = com_block_conf_quant_larger_trt(Z, Y, block,
    #                                          k.vec = (n - sum(Z) + 1) : n,
    #                                          methods.list.all, weight, opt.method,
    #                                          stat.null, null.max, tol, alpha = alpha / 2)

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

    # ci.control = com_block_conf_quant_larger_trt(Z, Y, block,
    #                                              k.vec = (n - sum(Z) + 1) : n,
    #                                              methods.list.all, weight, opt.method,
    #                                              stat.null, null.max, tol, alpha = alpha / 2)
    ci.all = sort( c(ci.treat, ci.control))
    return(ci.all)
  }
}






