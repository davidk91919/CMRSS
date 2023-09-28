
### Function for plotting heatmap, corresponding to multiple k(quantile of individual treatment effect) and c(upper bound on the null)

heatmap_pval = function(Z, Y, 
                        k.vec, ## vector corresponds to quantiles of interest
                        c.vec, ## vector corresponds to upper bounds on nulls
                        block,
                        methods.list.all,
                        weight = NULL,
                        stat.null = NULL,
                        null.max = 10^5,
                        Z.perm.all = NULL){
  require(CMRSS)
  require(pheatmap)
  
  l = length(k.vec)
  v = length(c.vec)
  
  res.mat = matrix(0, nrow = l, ncol = v)
  for(i in 1 : l){
    for(k in 1 : v)
      res.mat[i,k] = pval_comb(Z = Z, 
                               Y = Y, 
                               k = k.vec[i],
                               c = c.vec[k],
                               block = block,
                               methods.list.all = methods.list.all,
                               weight = weight, 
                               stat.null = stat.null,
                               null.max = null.max,
                               Z.perm.all = Z.perm.all,
                               statistic = FALSE)
  }  
  colnames(res.mat) = round(c.vec, 2)
  rownames(res.mat) = k.vec
  
  pheatmap(res.mat,
           cluster_rows = FALSE,
           cluster_cols = FALSE)    ## omitted clustering. can add an option for this
}



## function for power simulation; similar to 270p of the book. 
## working code; need to determine whether we are doing power simulation or try to get exact p-value under the model.
## for binary treatment, multiple(possibly not equal size of) strata
## now for simple additive model with tau vector, 
## and comparing difference in mean stat, wilcoxon stat and our new statistic

power_sim_block = function(Z, Y,
                           block,
                           k, c,
                           asym = TRUE, 
                           tau = seq(0, 1, len = 3), # treatment effect for alternative distribution; will be a vector, corresponding to different additive cases
                           methods.list.all,
                           weight = NULL, 
                           stat.null = NULL,
                           Z.perm.all = NULL,
                           null.max = 10^4, # length of the null distribution
                           iter.sim = 100, # number for simulations
                           alpha = 0.05){
  
  if(!is.factor(block)){
    block = as.factor(block)
  }  

  ## STEP 1: generating null distribution, and record threshold ##
  
  thres = rep(NA, 2)
  
  N = length(Z)
  B = length(levels(block))
  block.levels = levels(block)
  
  nb = rep(NA, B)
  mb = rep(NA, B)
  for(i in 1:B) {
    nb[i] = sum(block == block.levels[i])
    mb[i] = sum( Z[block == block.levels[[i]]] )
  } 
  indx = c(0, cumsum(nb) )
  
  new_block = as.factor( rep(1 : B, nb) )
  
  if(asym == TRUE){  
    fit_t = lm_robust(Y ~ Z, fixed_effects = ~block)
    thres[1] = confint(fit_t, level = (1 - alpha * 2) )[2]  #asymptotic t-stat
  } else {
    
    t.null.dist = rep(NA, null.max)
    
    ## generating potential outcomes; now, just for additive model; on the ** null distribution ** 
    
    Y0 = rep(NA, N)
    for(i in 1:B){
      trt.b = Z[block == block.levels[i]]
      Y.temp = Y[block == block.levels[i]]
      Y.temp[trt.b == 1] = Y.temp[trt.b == 1]       # null is fixed in this case(?)
      Y0[ (indx[i] + 1) : indx[i+1] ] = Y.temp
    }
    Y1 = Y0
    
    ## generating assignment matrix  
    
    if(is.null(Z.perm.all)){
      Z.perm = matrix(nrow = N, ncol = null.max)
      for(iter in 1 : null.max){
        for(i in 1:B){
          Z.perm[(indx[i] + 1 ):indx[i+1],iter] = sample( c(rep(1, mb[i]), rep(0, nb[i] - mb[i])) )
        }
      }
    }    else {       
      Z.perm = Z.perm.all
    }
    
    for(k in 1 : null.max){
      Z.iter = Z.perm[, k]
      Y.iter = Y0 * (1 - Z.iter) + Y1 * Z.iter
      t.null.dist[k] = lm_robust(Y.iter ~ Z.iter, fixed_effects = ~new_block)$statistic
    }
    
    thres[1] = sort(t.null.dist, decreasing = TRUE)[ floor(null.max * alpha) + 1]     
  }
  
  ## for new statistic
  
  ms_list = mu_sigma_list(Z, block = block, methods.list.all = methods.list.all)
  comb.null.dist = com_null_dist_block(Z, block = block, 
                                       methods.list.all = methods.list.all, 
                                       mu_sigma_list = ms_list,
                                       null.max = null.max)
  thres[2] = sort(comb.null.dist, decreasing = TRUE)[ floor(null.max * alpha) + 1]
  
  ## STEP 2: check the power ##
  ### this will come from the true distribution; true treatment effect we assumed ###
  
  reject.all = matrix(NA, nrow = 2, ncol = iter.sim)
  result = matrix(NA, nrow = length(tau), ncol = 2 )
  
  for(o in 1 : length(tau)){ 
    reject.all = matrix(NA, nrow = 2, ncol = iter.sim)
    
    for(iter in 1 : iter.sim){
      
      tau.sim = tau[o]
      temp = rep(NA, 2)
      
      #generate assignment(Z) corresponding to block structure
      
      Z.sim = rep(NA, N)
      for (i in 1 : length(nb) ){
        Z.sim[ (indx[i] + 1) : indx[i+1] ] = sample(rep(c(0,1), times = c( nb[i] - mb[i] , 
                                                                           mb[i]  ) ), 
                                                    nb[i], replace = FALSE)
      }
      Z_block = as.factor( rep(1 : B, nb) )
      block.levels.sim = levels(Z_block)
      
      # generate Y
      
      Y0.sim = rep(NA, N)
      for(i in 1:B){
        trt.b = Z.sim[Z_block == block.levels.sim[i]]
        Y.temp = Y[Z_block == block.levels.sim[i]]
        Y.temp[trt.b == 1] = Y.temp[trt.b == 1] - tau.sim   
        Y0.sim[ (indx[i] + 1) : indx[i+1] ] = Y.temp
      }
      Y1.sim = Y0.sim + tau.sim
      Y.sim = Z.sim * Y1.sim + (1 - Z.sim) * Y0.sim
      
      # calculate test statistics
      
      temp[1] = lm_robust(Y.sim ~ Z.sim, fixed_effects = ~Z_block)$statistic
      temp[2] = pval_comb(Z = Z.sim, 
                          Y = Y.sim,
                          k = k, 
                          c = c, 
                          block = Z_block,
                          methods.list.all = methods.list.all,
                          null.max = 10^3,
                          statistic = TRUE)[2]
      
      reject.all[,iter] = as.numeric( temp - thres > 0)
      print(iter)
    }
    result[o,] = rowMeans(reject.all)
  }  
  print(result)
}
