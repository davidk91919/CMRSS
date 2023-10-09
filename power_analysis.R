###############
### Function for plotting heatmap, corresponding to multiple k(quantile of individual treatment effect) and c(upper bound on the null)
###############

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

###############
## Function for power simulation, treating estimated quantile treatment effects as a true treatment effects. 
## rq function in quantreg package is used for estimating quantile treatment effect.
## This is for binary treatment, multiple(possibly not equal size of) strata
## usual difference in mean test statistic and our new test statistic are compared
###############


# estimate quantile effect by quantreg package and treat this as a true effect

library(quantreg)

Y.obs = pbs1$R01TMCRET
Z.obs = pbs1$trt
block.obs = factor(pbs1$STRA_BLOCK)

N = length(Z.obs)
B = length(levels(block.obs))
block.levels = levels(block.obs)

nb = rep(NA, B)
mb = rep(NA, B)
for(i in 1:B) {
  nb[i] = sum(block.obs == block.levels[i])
  mb[i] = sum( Z.obs[block.obs == block.levels[[i]]] )
} 
indx = c(0, cumsum(nb) )

Z_block = as.factor( rep(1 : B, nb) )
block.levels = levels(Z_block)

Y0_tmp = rep(NA, N)
Y1_tmp = rep(NA, N)
trt_eff = rep(NA, N)

for(u in 1 : B){
  trt.b = Z.obs[Z_block == block.levels[u]]
  Y.b = Y.obs[Z_block == block.levels[u]]
  prob.seq.tmp = seq(0,1, length.out = length(Y.b))
  tmp.data = data.frame(y = Y.b, trt = trt.b)
 
  m.temp = rq(y ~ trt, tau = prob.seq.tmp, data = tmp.data)
  trt_eff[ (indx[u] + 1) : indx[u+1] ] = m.temp$coef[2,]
  Y0_tmp[ (indx[u] + 1) : indx[u+1] ] = Y.b - trt.b * m.temp$coef[2,]
  Y1_tmp[ (indx[u] + 1) : indx[u+1] ] = Y.b + (1 - trt.b) * m.temp$coef[2,]
} # note Y0 and Y1 is also arranged here by block orders

any(is.na(Y0_tmp))
any(is.na(Y1_tmp))
min(trt_eff) > 0 # note there exists some negative treatment effect under 75%
quantile(trt_eff)

# function for power simulation, treating estimated quantile treatment effect as true

power_sim_block_quant = function(Z, 
                                 Y0, Y1, # potential outcome vectors; will estimate by ci.qte object
                                 block,
                                 k, 
                                 c,
                                 methods.list.all,
                                 weight = NULL, 
                                 iter.sim = 100, # number for simulations
                                 alpha = 0.05
){
  
  if(!is.factor(block)){
    block = as.factor(block)
  }  
  
  N = length(Z)
  B = length(levels(block))
  #  block.levels = levels(block)
  block.levels = levels(block)
  
  nb = rep(NA, B)
  mb = rep(NA, B)
  for(i in 1:B) {
    nb[i] = sum(block == block.levels[i])
    mb[i] = sum( Z[block == block.levels[[i]]] )
  } 
  indx = c(0, cumsum(nb) )
  Z_block = as.factor( rep(1 : B, nb) )
  block.levels = levels(Z_block)
  
  result = rep(NA, 2)
  
  reject.all = matrix(NA, nrow = 2, ncol = iter.sim)
  
  for(iter in 1 : iter.sim){
    temp = rep(NA, 2)
    
    #generate assignment(Z) corresponding to block structure
    
    Z.sim = rep(NA, N)
    for (i in 1 : length(nb) ){
      Z.sim[ (indx[i] + 1) : indx[i+1] ] = sample(rep(c(0,1), times = c( nb[i] - mb[i] , 
                                                                         mb[i]  ) ), 
                                                  nb[i], replace = FALSE)
    }
    
    Y.sim = Z.sim * Y1 + (1 - Z.sim) * Y0

    # calculate test statistics and compare their p-values
    
    temp[1] = lm_robust(Y.sim ~ Z.sim, fixed_effects = ~Z_block)$p.value  
    temp[2] = pval_comb(Z = Z.sim, 
                        Y = Y.sim,
                        k = k, 
                        c = c, 
                        block = Z_block,
                        methods.list.all = methods.list.all,
                        null.max = 10^3,
                        statistic = FALSE)
    
    reject.all[,iter] = as.numeric( temp < alpha )
    print(iter)
  }
  result[1,] = rowMeans(reject.all)
  print(result)
}

# test code

power_sim_block_quant(Z = pbs1$trt, 
                      Y0 = Y0_tmp, 
                      Y1 = Y1_tmp,
                      block = factor(pbs1$STRA_BLOCK), 
                      k = length(pbs1$R01TMCRET),
                      c = 0,
                      methods.list.all = methods.list,
                      iter.sim = 100,
                      alpha = 0.05)
