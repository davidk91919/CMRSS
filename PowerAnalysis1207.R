library(QIoT)
library(RIQITE)
library(CMRSS)

## This is for 'dumb' simulation for both CRE(RIQITE) and SERE(QIoT), comparing to CMRSS

## data generation: replicating block and treatment structure of LC data

N = 5768
B = 44
nb =  c(27,   16,   31,   51,   38,   31,   18,   18,   55,   28,   52,   31,   36,
        43,   28,   36,   32,   23,   38,   26,   63,   31,   11,   25,   14,   29,
        24,   62,  176,  380,  362,  506,  122,  371,  192,  398, 1273,   82,  288,
        266,  227,   89,   58,   61)
mb = c(15,   7,  15,  30,  23,  19,   9,  10,  31,  17,  31,  19,  22,  26,  18,  21,
       19,  14,  22,  16,  36,  18,   6,  15,   9,  18,  14,  37,  89, 188, 182, 252,
       75, 221, 115, 239, 761,  55, 191, 178, 152,  54,  38,  41)
indx = c(0, cumsum(nb) )


## CRE

# listing methods for both RIQITE and CMRSS

m = sum(mb)

m.vec = floor(seq(2, m, len = 5))
alpha = 0.05
k = floor(N * 0.95)
c = 0

method.list.1 = list(name = "Stephenson", scale = FALSE,
                     s = m.vec[1])
method.lists.1 = list(list(method.list.1))

method.list.2 = list(name = "Stephenson", scale = FALSE,
                     s = m.vec[2])
method.lists.2 = list(list(method.list.2))

method.list.3 = list(name = "Stephenson", scale = FALSE,
                     s = m.vec[3])
method.lists.3 = list(list(method.list.3))

method.list.4 = list(name = "Stephenson", scale = FALSE,
                     s = m.vec[4])
method.lists.4 = list(list(method.list.4))

method.list.5 = list(name = "Stephenson", scale = FALSE,
                     s = m.vec[5])
method.lists.5 = list(list(method.list.5))

H = 5
s = 1

method.list.all.1 <- list()
method.list.all.2 <- list()
method.list.all.3 <- list()
method.list.all.4 <- list()
method.list.all.5 <- list()

methods.list = list(method.lists.1,
                    method.lists.2,
                    method.lists.3,
                    method.lists.4,
                    method.lists.5)


for (j in 1 : H){
  for(i in 1 : s){
    methods.list[[j]][[i]] = list(name = "Polynomial",
                                  r = m.vec[j],
                                  std = TRUE,
                                  scale = TRUE)
  }
}

# function specification



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
  
  #  N = length(Z)
  N = length(Y0)
  B = length(levels(block))
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
    
    Y.tmp = Z.sim * Y1 + (1 - Z.sim) * Y0
    Y.sim = sample(Y.tmp, size = N, replace = FALSE)   ## added for permute outcome indices before analysis
    
    
    # calculate test statistics and compare their p-values
    
    temp[1] = pval_quantile(Z = Z.sim,
                            Y = Y.sim,
                            k = k,
                            c = c,
                            method.list = list(name = "Stephenson", s = 150),
                            nperm = 10^3)
    
    temp[2] = pval_comb(Z = Z.sim,
                        Y = Y.sim,
                        k = k,
                        c = c,
                        block = Z_block,
                        methods.list.all = methods.list.all,
                        null.max = 10^3,
                        statistic = FALSE,
                        opt.method = "ILP_gurobi")
    
    reject.all[,iter] = as.numeric( temp < alpha )
    print(iter)
  }
  result = rowMeans(reject.all)
  print(result)
}

# generating 'dumb' potential outcome

Y0_tmp = rep(NA, N) 
Y1_tmp = rep(NA, N)
Y = rep(0, N)
#Y = runif(N) ## case for non-tied outcomes
dumb.prop = 0.1 ## proportion of effected units

eff.len = floor(N * dumb.prop)    ### since k on tau_k is about the whole unit
eff.size = rep(1, eff.len)
#eff.size = runif(n = eff.len, min = 0.5, max = 1) ### case for non-tied effects
trt.b = sample(rep( c(0,1), c(N - m, m),
                    size = N), replace = FALSE)

trt.idx = which(trt.b == 1)
neff.idx = sample(trt.idx, size = m - eff.len, replace = FALSE)
eff.idx = setdiff(trt.idx, neff.idx)

eff.b = trt.b
eff.b[-trt.idx] = 0              
eff.b[neff.idx] = 0
eff.b[eff.idx] = eff.size

Z.obs = trt.b
Y0_tmp = Y - trt.b * eff.b
Y1_tmp = Y + (1 - trt.b) * eff.b
quantile(Y1_tmp - Y0_tmp, probs = seq(0.3, 1, len = 23)) # check quantile treatment effect generated as desired

# implement power analysis on dumb simulation data


Z_block = factor(rep(1, N))

power_sim_block_quant(Z = Z.obs,
                      Y0 = Y0_tmp,
                      Y1 = Y1_tmp,
                      block = Z_block,
                      k = floor(0.95 * length(Y0_tmp) ),
                      #c = -2e-08,
                      c = 0,
                      methods.list.all = methods.list,
                      iter.sim = 5,
                      alpha = 0.05)



## SERE

Z_block = as.factor( rep(1 : B, nb) )
block.levels = levels(Z_block)

Y0_tmp = rep(NA, N)
Y1_tmp = rep(NA, N)
Z.obs = rep(NA, N)
dumb.prop = 0.1 ## hypothetical proportion of effect
eff.size = 1

for(u in 1 : B){
  Y.b = rep(1, nb[u])
  #  Y.b = runif(nb[u]) ## for no-tied outcomes
  eff.len = floor(nb[u] * dumb.prop)    ### since k on tau_k is about the whole unit
  
  trt.b = sample(rep( c(0,1), c(nb[u]- mb[u], mb[u]),
                      size = nb[u]), replace = FALSE)
  
  trt.idx = which(trt.b == 1)
  neff.idx = sample(trt.idx, size = mb[u] - eff.len, replace = FALSE)
  eff.idx = setdiff(trt.idx, neff.idx)
  
  eff.b = trt.b
  eff.b[-trt.idx] = 0
  eff.b[neff.idx] = 0
  eff.b[eff.idx] = rep(eff.size, len = length(eff.idx))
  #  eff.b[eff.idx] = seq(0.5, 1, len = length(eff.idx) ) # for no-tied effect
  #  eff.b = c( rep(0, length(trt.b) - eff.len), seq(0.0001, 1, len = eff.len))
    Z.obs[ (indx[u] + 1) : indx[u+1] ] = trt.b
    Y0_tmp[ (indx[u] + 1) : indx[u+1] ] = Y.b - trt.b * eff.b
    Y1_tmp[ (indx[u] + 1) : indx[u+1] ] = Y.b + (1 - trt.b) * eff.b
  
#  Z.obs[ lc_dat$st_site_block == levels(block.obs)[u]] = trt.b
#  Y0_tmp[ lc_dat$st_site_block == levels(block.obs)[u]] = Y.b
#  Y1_tmp[ lc_dat$st_site_block == levels(block.obs)[u]] = Y.b + eff.b #### if we know the original block structure
  
  
} # note Y0 and Y1 is also arranged here by block orders

any(is.na(Y1_tmp))
quantile(Y1_tmp - Y0_tmp, probs = seq(0.3, 1, len = 23))
quantile(Y1_tmp, probs = seq(0, 1, len = 23))
unique(eff.b)
length(unique(Y1_tmp))
length(unique(Y0_tmp))

# implement power analysis

power_sim_block_quant(Z = Z.obs,
                      Y0 = Y0_tmp,
                      Y1 = Y1_tmp,
                      block = Z_block,
                      k = floor(0.95 * length(Y0_tmp) ),
                      c = 0,
                      #c = -2e-08,
                      methods.list.all = methods.list,
                      iter.sim = 5,
                      alpha = 0.05)






