
library(QIoT)
library(CMRSS)
library(quantreg)

# methods
## for QIoT

s = 44

methods.list.qiot.1 = list()
methods.list.qiot.2 = list()
methods.list.qiot.3 = list()

for(i in 1 : s) {
  methods.list.qiot.1[[i]] <- list(name = "Stephenson",
                                 s = 6)
}

for(i in 1 : s) {
  methods.list.qiot.2[[i]] <- list(name = "Stephenson",
                                   s = 24)
}

for(i in 1 : s) {
  methods.list.qiot.3[[i]] <- list(name = "Stephenson",
                                   s = 55)
}

## for CMRSS

H = 10

m.seq = c(2,  87, 171, 255, 340, 424, 508, 593, 677, 761)

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

methods.list      <-     list(method.list.all.1, method.list.all.2,
                              method.list.all.3, method.list.all.4,
                              method.list.all.5, method.list.all.6,
                              method.list.all.7, method.list.all.8,
                              method.list.all.9, method.list.all.0)


for (j in 1 : H){
  for(i in 1 : s){
    methods.list[[j]][[i]] = list(name = "Polynomial",
                                  r = m.seq[j],
                                  std = TRUE,
                                  scale = TRUE)
  }
}




# power simulation function

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

  result = rep(NA, 4)

  reject.all = matrix(NA, nrow = 4, ncol = iter.sim)

  for(iter in 1 : iter.sim){
    temp = rep(NA, 4)

    #generate assignment(Z) corresponding to block structure

    Z.sim = rep(NA, N)
    for (i in 1 : length(nb) ){
      Z.sim[ (indx[i] + 1) : indx[i+1] ] = sample(rep(c(0,1), times = c( nb[i] - mb[i] ,
                                                                         mb[i]  ) ),
                                                  nb[i], replace = FALSE)
    }

    Y.sim = Z.sim * Y1 + (1 - Z.sim) * Y0

    # calculate test statistics and compare their p-values

    temp[1] = pval_quantile_scre(Z = Z.sim,
                                 Y = Y.sim,
                                 k = k,
                                 c = c,
                                 block = Z_block,
                                 method.list.all = methods.list.qiot.1,
                                 opt.method = "ILP_gurobi")$lower
    temp[2] = pval_quantile_scre(Z = Z.sim,
                                 Y = Y.sim,
                                 k = k,
                                 c = c,
                                 block = Z_block,
                                 method.list.all = methods.list.qiot.2,
                                 opt.method = "ILP_gurobi")$lower
    temp[3] = pval_quantile_scre(Z = Z.sim,
                                 Y = Y.sim,
                                 k = k,
                                 c = c,
                                 block = Z_block,
                                 method.list.all = methods.list.qiot.3,
                                 opt.method = "ILP_gurobi")$lower
    temp[4] = pval_comb(Z = Z.sim,
                        Y = Y.sim,
                        k = k,
                        c = c,
                        block = Z_block,
                        methods.list.all = methods.list.all,
                        null.max = 10^4,
                        statistic = FALSE)

    reject.all[,iter] = as.numeric( temp < alpha )
    print(iter)
  }
  result = rowMeans(reject.all)
  print(result)
}



## some 'dumb' simulation
# with the same block structure


N = 5768
B = 44
nb =  c(27,   16,   31,   51,   38,   31,   18,   18,   55,   28,   52,   31,   36,
43,   28,   36,   32,   23,   38,   26,   63,   31,   11,   25,   14,   29,
24,   62,  176,  380,  362,  506,  122,  371,  192,  398, 1273,   82,  288,
266,  227,   89,   58,   61)
mb = c(15,   7,  15,  30,  23,  19,   9,  10,  31,  17,  31,  19,  22,  26,  18,  21,
19,  14,  22,  16,  36,  18,   6,  15,   9,  18,  14,  37,  89, 188, 182, 252,
75, 221, 115, 239, 761,  55, 191, 178, 152,  54,  38,  41)

Z_block = as.factor( rep(1 : B, nb) )
block.levels = levels(Z_block)

Y0_tmp = rep(NA, N)
Y1_tmp = rep(NA, N)
trt_eff = rep(NA, N)
for(u in 1 : B){
  trt.b = Z.obs[Z_block == block.levels[u]]
  Y.b = rep(0, length(trt.b))
  eff.len = floor(length(trt.b) * 0.1)
  eff.b = rep(c(0, 1), c( length(trt.b) - eff.len , eff.len))

  Y0_tmp[ (indx[u] + 1) : indx[u+1] ] = Y.b
  Y1_tmp[ (indx[u] + 1) : indx[u+1] ] = Y.b + eff.b
} # note Y0 and Y1 is also arranged here by block orders

quantile(Y1_tmp - Y0_tmp, probs = seq(0,1, len = 23))

power_sim_block_quant(Z = Z.obs,
                      Y0 = Y0_tmp,
                      Y1 = Y1_tmp,
                      block = Z_block,
                      k = floor(0.95 * length(Y0_tmp) ),
                      c = 0,
                      methods.list.all = methods.list,
                      iter.sim = 10,
                      alpha = 0.05)










##########
## original code for power analysis
##########



Y.obs = lc_dat$R01ATMCRET
Z.obs = lc_dat$trt
block.obs = factor(lc_dat$st_site_block)

N = length(Z.obs)
B = s = length(levels(block.obs))
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
quantile(trt_eff,
         probs = c(0, 0.1, 0.3, 0.5, 0.8, 0.9, 0.93, 0.95, 0.97, 1))

power_sim_block_quant(Z = Z.obs,
                      Y0 = Y0_tmp,
                      Y1 = Y1_tmp,
                      block = Z_block,
                      k = floor(length(Y0_tmp) ),
                      c = 0,
                      methods.list.all = methods.list,
                      iter.sim = 6,
                      alpha = 0.05)


