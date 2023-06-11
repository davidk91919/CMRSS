#####################################################
# This is a script of simulations for both CRE and SERE
#####################################################

################################
# Required Source and Libraries
################################
#source("Fun_multiple_stephenson_0126.R")
library(CMRSS)
library(gurobi)
library(psych)
library(randomizr)
library(Matrix)
library(skewt)

##################################
# Simulation Codes for CRE(updated; for both polynomial and Stephenson, also for combined in both ways)
##################################


n = 150         # total number of units
m = 0.5 * n     # number of treated units
nperm = 10^4    # total number of permutation in simulation
alpha = 0.1     
k = ceiling(n * 0.9)    ## testing for 90% quantile treatment effect
c = 0


s.vec = c(2, 4, 10, 20, 30, 60, 90, 120, 150, 200, 250)    ## parameter vector for both Stephenson and Polynomial rank score statistic

H = length(s.vec)
stat.null.all <- matrix(NA, nrow = nperm, ncol =  ( 2* length(s.vec) + 2 ) )
Z.perm = assign_CRE(n, m, nperm)
thres.all <- rep(NA, (2 * length(s.vec) + 2) )    ## generate threshold vector for simulated test statistics

for(i in 1: H ){
  stat.null.all[, i] = null_dist(n, m, 
                                 method.list = list(name = "Stephenson", 
                                                    s = s.vec[i], 
                                                    scale = TRUE), 
                                 nperm = 10^5, Z.perm = Z.perm)
  # find threshold such that p-value <= alpha <===> test stat > threshold #
  thres.all[i] = thres = sort(stat.null.all[, i], decreasing = TRUE)[ floor(nperm * alpha) + 1 ]
  
  stat.null.all[, H + i] = null_dist(n, m,
                                     method.list = list(name = "Polynomial",
                                                        r = s.vec[i],
                                                        std = TRUE,
                                                        scale = TRUE),
                                     nperm = 10^5, Z.perm = Z.perm)
  thres.all[H + i] = sort(stat.null.all[, H + i], decreasing = TRUE)[ floor(nperm * alpha) + 1 ]
}

## threshold for combinded Stephenson
stat.null.all[, 2 * length(s.vec) + 1] = apply( stat.null.all[, 1 : (length(s.vec) - 2)  ], 1, max)
thres.all[2 * length(s.vec) + 1] = thres = sort(stat.null.all[, 2 * length(s.vec) + 1], decreasing = TRUE)[ floor(nperm * alpha) + 1]

# threshold for combinded Polynomial
stat.null.all[, 2 * length(s.vec) + 2] = apply( stat.null.all[, (H + 1) : (H + length(s.vec) ) ], 1, max)
thres.all[2 * length(s.vec) + 2] = thres = sort(stat.null.all[, 2 * length(s.vec) + 2], decreasing = TRUE)[ floor(nperm * alpha) + 1]


iter.max = 1000
sig = seq(0,5, by = 0.5)    ## sequence of sigma; variance of simulated outcome
r.sim = matrix(NA, nrow = length(sig), ncol = (2 * length(s.vec) + 2) ) ## matrix resulted by simulaton. each elements corresponds to each cases - sigma and method
eff.r.sim = rep(NA, length(sig))      ## true effect vector corresponding to each sigma
p = n - k 

for(j in 1 : length(sig)){
  reject.all = matrix(NA, nrow = (2 * length(s.vec) + 2), ncol = iter.max)
  for(iter in 1:iter.max){
    
    ## data generating model
    #   Y1 = rnorm(n) + 2 * ( runif(n) <= 0.6 )       #a1
    #   Y0 = rnorm(n)
    
    #   Y1 = 0.8 * rnorm(n) + 1 * ( runif(n) <= 0.8 )
    #   Y0 = rnorm(n)
    #   Y0 = rnorm(n)
    #   Y1 = Y0 + sig[j] * rnorm(n)
    #   Y1 = sig[j] * rnorm(n)
    
    
      Y1 = 0.8 * rnorm(n) + 2 * ( runif(n) <= 0.6 )
      Y0 = rnorm(n)
    
    Z = assign_CRE(n, m, 1)
    Y = Z * Y1 + (1-Z) * Y0
    rbind(Z, Y)
    
    ind.sort.treat = sort_treat(Y, Z)
    
    # for a single s
    min.stat.all = rep(NA, (2 * length(s.vec)) )
    
    for(u in 1 : length(s.vec) ){
      min.stat.all[u] = min_stat(Z, Y, k, c, 
                                 method.list = list(name = "Stephenson", 
                                                    s = s.vec[u],
                                                    scale = TRUE),
                                 ind.sort.treat = ind.sort.treat)
      min.stat.all[u + H] = min_stat(Z, Y, k, c, 
                                     method.list = list(name = "Polynomial",
                                                        r = s.vec[u],
                                                        std = TRUE,
                                                        scale = TRUE), 
                                     ind.sort.treat = ind.sort.treat)
    }
    
    min.stat.all[2 * length(s.vec) + 1] <- max(min.stat.all[1 : ( length(s.vec) - 2 )])
    min.stat.all[2 * length(s.vec) + 2] <- max(min.stat.all[(H + 1) : (H + length(s.vec))])
    reject.all[, iter] = as.numeric( min.stat.all - thres.all > 0 )
    
    print(iter)  
  }
  
  rowMeans(reject.all)
  r.sim[j, ] = rowMeans(reject.all)
  eff.r.sim[j] = sort(Y1 - Y0, decreasing = TRUE)[p]
}

r.sim
eff.r.sim




##################################
# Simulation Codes for SERE
##################################
##################################
# for equal size of strata (n = 50 for each strata)
##################################
## Simulation Data and listing methods

s = 100     ## number of strata
n = 50      ## number of units in each strata
m = 0.5 * n ## number of treated units
N = s * n   
k = ceiling(0.95 * N) ## k corresponding to the null H_0: \tau_{k} \leq c
p = N - k
c = 0

simulation_Z <- as.data.frame(block.random(n = N, c(block = s))) ## only needed for equal strata

simulation_Z$block <- as.factor(simulation_Z$block)
Z = block_ra(blocks = simulation_Z$block)
Z_block <- simulation_Z$block

### listing list of methods

method.list.all.1 <- list()
method.list.all.2 <- list()
method.list.all.3 <- list()
method.list.all.4 <- list()
method.list.all.5 <- list()
method.list.all.6 <- list()
method.list.all.7 <- list()
method.list.all.8 <- list()
method.list.all.9 <- list()
method.list.all.0 <- list()


comb.methods.list.all <- list(method.list.all.1, method.list.all.2,
                         method.list.all.3, method.list.all.4,
                         method.list.all.5, method.list.all.6,
                         method.list.all.7, method.list.all.8,
                         method.list.all.9, method.list.all.0)

prev.methods.list.all <- list(method.list.all.1, method.list.all.2,
                         method.list.all.3, method.list.all.4,
                         method.list.all.5, method.list.all.6,
                         method.list.all.7, method.list.all.8,
                         method.list.all.9, method.list.all.0)

H = h = length(comb.methods.list.all)
r.seq = ceiling(seq(2, n, len = H))
# r = c(2, 8, 13, 18, 24, 29, 34, 40, 45, 50) here 

## specifying methods for Polynomial rank score statistics

for (j in 1 : H){
  for(i in 1 : s){
    comb.methods.list.all[[j]][[i]] = list(name = "Polynomial",
                                           r = r.seq[j],
                                           std = TRUE,
                                           scale = TRUE)
  }
}

## specifying methods for Stephenson(Wilcoxon) rank score statistics

for(i in 1 : s){
  prev.methods.list.all[[1]][[i]] = list(name = "Wilcoxon",
                                         scale = TRUE)
}


for (j in 2 : H){
  for(i in 1 : s){
    prev.methods.list.all[[j]][[i]] = list(name = "Stephenson",
                                           s = r.seq[j],
                                           scale = TRUE)
  }
}


## generating alpha-level thresholds for simulation

nperm = 10^5  ## decreased, since R session is getting down for 10^5
thres.all = rep(0, 2 * h + 1)
alpha = 0.05 # threshold
ms_list = mu_sigma_list(Z, block = Z_block, methods.list.all = comb.methods.list.all)

comb.null.dist = com_null_dist_block(Z, Z_block, comb.methods.list.all, mu_sigma_list = ms_list,
                                     null.max = nperm)
thres.all[2 * h + 1] = aa = sort(comb.null.dist, decreasing = TRUE)[ floor(nperm * alpha) + 1]

for(u in 1 : h){
  null.dist.1 = null_dist_block(Z, Z_block, prev.methods.list.all[[u]],
                                null_max = nperm)
  thres.all[u] = sort(null.dist.1, decreasing = TRUE)[ floor(nperm * alpha) + 1]
}
for(j in 1 : h){
  methods.list.all = list(comb.methods.list.all[[j]])
  ms_list_tmp = list(mean = ms_list$mean[j],
                     sigma = ms_list$sigma[j])
  null.dist.2 = com_null_dist_block(Z, Z_block, methods.list.all, mu_sigma_list = ms_list_tmp,
                                    null.max = nperm)
  thres.all[h + j] = sort(null.dist.2, decreasing = TRUE)[ floor(nperm * alpha) + 1]
}

## check whether there exists some overwritten variable names
#k
#p
#N
#n


## simulation implement
# null in here ; H_0: \tau_{(900)} \leq 0 

iter.max = 10^2 ## total number for simulations
sig = seq(0,5, by = 0.5)
p.sim.eq.95.loc = matrix(NA, nrow = length(sig), ncol = 2 * H + 1)   ## matrix resulted by simulaton. each elements corresponds to each cases - sigma and method
eff.p.sim.eq.95.loc = rep(NA, length(sig))  ## vector for true effect size for each cases(each sigma)

for(j in 1 : length(sig)){
  reject.all = matrix(0, nrow = 2 * H + 1, ncol = iter.max)
  for(iter in 1 : iter.max){
    
    #  temp = rep(0, h + 1)
    temp = rep(0, 2 * h + 1)
    
    ## generate test assignment / outcome
    simulation_Z <- as.data.frame(block.random(n = N, c(block = s)))
    simulation_Z$block <- as.factor(simulation_Z$block)
    Z = block_ra(blocks = simulation_Z$block)
    Z_block <- simulation_Z$block
    #Y1 = rnorm(N) - 1
    #Y0 = rnorm(N)
    
    #  Y1 = rnorm(N) + 4 * ( runif(N) <= 0.5 )
    #  Y0 = rnorm(N)  
    
    #  Y1 = rnorm(N, mean = 0.5, sd = 20) 
    #  Y0 = rnorm(N)
    
    Y0 = rep(NA, N)
    Y1 = rep(NA, N)
    for(i in 1 : s){
      #    n = length(Z[Z_block == i])
      #    tmp = rnorm(n[i])
      Y0[Z_block == i] = rnorm(n)
      Y1[Z_block == i] = rnorm(n) * sig[j] + 0.5
    }
    Y = Z * Y1 + (1-Z) * Y0
    
    ## calculate min_stat for each cases
    
    coeflists.comb <- comb_matrix_block(Z, Y, block = Z_block,
                                        c = c, 
                                        methods.list.all = comb.methods.list.all)
    
    temp[2 * h+1] = Gurobi_sol_com(Z, block = Z_block, weight = NULL, 
                                   coeflists.comb, p, mu_sigma_list = ms_list, exact = TRUE)$obj
    for (i in 1:h){
      coeflists      <- comb_matrix_block(Z, Y, block = Z_block,
                                          c = c, 
                                          methods.list.all = list(comb.methods.list.all[[i]])
      )
      ms_list.tmp = list(mean = ms_list$mean[i], sigma = ms_list$sigma[i])
      
      temp[i] = QIoT::min_stat_block(Z, Y, block = Z_block, k, c, 
                                     method.list.all = prev.methods.list.all[[i]], 
                                     opt.method = "ILP_gurobi" )$lower
      temp[i + h] = Gurobi_sol_com(Z, block = Z_block, weight = NULL, 
                                   coeflists, p, mu_sigma_list = ms_list.tmp, exact = TRUE)$obj
    }
    
    reject.all[,iter] = as.numeric( temp - thres.all > 0)
    
    print(iter)
  }
  p.sim.eq.95.loc[j,] = rowMeans(reject.all)
  eff.p.sim.eq.95.loc[j] = sort(Y1 - Y0, decreasing = TRUE)[p]
}

##################################
# for equal size of strata (n = 50 for each strata)
##################################
## Simulation Data and listing methods

s = 100
n = sample(20 : 60, size = s, replace = TRUE)   ## units for each strata; now this is unequal
n.temp = c(0, cumsum(n))
N = sum(n)
Z = rep(0, N)
for (i in 1 : (length(n.temp) - 1) ){
  Z[(n.temp[i] + 1) : n.temp[i + 1] ] = sample(rep(c(0,1), times = c(ceiling(n[i] / 2 ), n[i] - ceiling(n[i] / 2)   ) ), n[i], replace = FALSE)
}
Z_block = as.factor( rep(1:s, times = n) )
k = ceiling(0.95 * N)
p = N - k
c = 0

# specifying methods

method.list.all.1 <- list()
method.list.all.2 <- list()
method.list.all.3 <- list()
method.list.all.4 <- list()
method.list.all.5 <- list()
method.list.all.6 <- list()
method.list.all.7 <- list()
method.list.all.8 <- list()
method.list.all.9 <- list()
method.list.all.10 <- list()

comb.methods.list.all <- list(method.list.all.1, method.list.all.2,
                              method.list.all.3, method.list.all.4,
                              method.list.all.5, method.list.all.6,
                              method.list.all.7, method.list.all.8,
                              method.list.all.9, method.list.all.10)
prev.methods.list.all <- list(method.list.all.1, method.list.all.2,
                              method.list.all.3, method.list.all.4,
                              method.list.all.5, method.list.all.6,
                              method.list.all.7, method.list.all.8,
                              method.list.all.9, method.list.all.10)

H = h = length(comb.methods.list.all)

## if we want to use different hyperparameters('r' or 's') in each strata for same method

r.seq = list()
r.seq = ceiling(seq(2, 40, len = H)) 

for(j in 1 : H){
  for(i in 1 : s){
    comb.methods.list.all[[j]][[i]] = list(name = "Polynomial",
                                           r = r.seq[j],
                                           std = TRUE,
                                           scale = FALSE)
  }
}


for(i in 1 : s){
  prev.methods.list.all[[1]][[i]] = list(name = "Wilcoxon",
                                         scale = TRUE)
}


for (j in 2 : H){
  for(i in 1 : s){
    prev.methods.list.all[[j]][[i]] = list(name = "Stephenson",
                                           s = r.seq[j],
                                           scale = TRUE)
  }
}

## if we want to use ** same ** hyperparameters ('r' or 's') in each strata for same method
### this will be mainly used

max.n = max(n)
r.seq = ceiling(seq(2, max.n, length.out = H))

for (j in 1 : H){
  for(i in 1 : s){
    comb.methods.list.all[[j]][[i]] = list(name = "Polynomial",
                                           r = r.seq[j],
                                           std = TRUE,
                                           scale = FALSE)
  }
}


for(i in 1 : s){
    prev.methods.list.all[[1]][[i]] = list(name = "Wilcoxon",
                                           scale = TRUE)
}


for (j in 2 : H){
  for(i in 1 : s){
    prev.methods.list.all[[j]][[i]] = list(name = "Stephenson",
                                           s = r.seq[j],
                                           scale = TRUE)
  }
}


# calculating threshold values

nperm = 10^5  
thres.all = rep(0, 2 * h + 1)
alpha = 0.05 # threshold
ms_list = mu_sigma_list(Z, block = Z_block, methods.list.all = comb.methods.list.all)

comb.null.dist = com_null_dist_block(Z, Z_block, comb.methods.list.all, mu_sigma_list = ms_list,
                                     null.max = nperm)
thres.all[2 * h + 1] = aa = sort(comb.null.dist, decreasing = TRUE)[ floor(nperm * alpha) + 1]

for(u in 1 : h){
  null.dist.1 = null_dist_block(Z, Z_block, prev.methods.list.all[[u]],
                                null_max = nperm)
  thres.all[u] = sort(null.dist.1, decreasing = TRUE)[ floor(nperm * alpha) + 1]
}
for(j in 1 : h){
  methods.list.all = list(comb.methods.list.all[[j]])
  ms_list_tmp = list(mean = ms_list$mean[j],
                     sigma = ms_list$sigma[j])
  null.dist.2 = com_null_dist_block(Z, Z_block, methods.list.all, mu_sigma_list = ms_list_tmp,
                                    null.max = nperm)
  thres.all[h + j] = sort(null.dist.2, decreasing = TRUE)[ floor(nperm * alpha) + 1]
}

# simulation implement

iter.max = 10^2
sig = seq(0,5, by = 0.5)
p.sim.uneq.95.loc = matrix(NA, nrow = length(sig), ncol = 2 * H + 1)    ## matrix resulted by simulaton. each elements corresponds to each cases - sigma and method
eff.p.sim.uneq.95.loc = rep(NA, length(sig))  ## vector for true effect size. 

for(j in 1 : length(sig)){
  reject.all = matrix(0, nrow = 2 * H + 1, ncol = iter.max)
  for(iter in 1 : iter.max){
    
    #  temp = rep(0, h + 1)
    temp = rep(0, 2 * h + 1)
    
    ## generate test assignment / outcome
    Z = rep(NA, N)
    for (i in 1 : length(n) ){
      Z[(n.temp[i] + 1) : n.temp[i + 1] ] = sample(rep(c(0,1), times = c(ceiling(n[i] / 2 ), n[i] - ceiling(n[i] / 2)   ) ), n[i], replace = FALSE)
    }
    Z_block = as.factor( rep(1:s, times = n) )
    #Y1 = rnorm(N) - 1
    #Y0 = rnorm(N)
    
    #  Y1 = rnorm(N) + 4 * ( runif(N) <= 0.5 )
    #  Y0 = rnorm(N)  
    
    #  Y1 = rnorm(N, mean = 0.5, sd = 20) 
    #  Y0 = rnorm(N)
    
    Y0 = rep(NA, N)
    Y1 = rep(NA, N)
    for(i in 1 : s){
      #    n = length(Z[Z_block == i])
      #    tmp = rnorm(n[i])
      Y0[Z_block == i] = rnorm(n[i])
      Y1[Z_block == i] = rnorm(n[i]) * sig[j] + 0.5
    }
    Y = Z * Y1 + (1-Z) * Y0
    
    ## calculate min_stat for each cases
    
    coeflists.comb <- comb_matrix_block(Z, Y, block = Z_block,
                                        c = c, 
                                        methods.list.all = comb.methods.list.all)
    
    temp[2 * h+1] = Gurobi_sol_com(Z, block = Z_block, weight = NULL, 
                                   coeflists.comb, p, mu_sigma_list = ms_list, exact = TRUE)$obj
    for (i in 1:h){
      coeflists      <- comb_matrix_block(Z, Y, block = Z_block,
                                          c = c, 
                                          methods.list.all = list(comb.methods.list.all[[i]])
      )
      ms_list.tmp = list(mean = ms_list$mean[i], sigma = ms_list$sigma[i])
      
      temp[i] = QIoT::min_stat_block(Z, Y, block = Z_block, k, c, 
                                     method.list.all = prev.methods.list.all[[i]], 
                                     opt.method = "ILP_gurobi" )$lower
      temp[i + h] = Gurobi_sol_com(Z, block = Z_block, weight = NULL, 
                                   coeflists, p, mu_sigma_list = ms_list.tmp, exact = TRUE)$obj
    }
    
    reject.all[,iter] = as.numeric( temp - thres.all > 0)
    
    print(iter)
  }
  p.sim.uneq.95.loc[j,] = rowMeans(reject.all)
  eff.p.sim.uneq.95.loc[j] = sort(Y1 - Y0, decreasing = TRUE)[p]
}




