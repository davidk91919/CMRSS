##### this is a building function for 'testing twice' for CMRSS.
##### this will finally generate less conservative p-value than before, by comparing all possible pairs of (n1, n0) = (a, n - k - 1 - a). 

switch_pval_comb <- function(Z, Y, k, c,
                             block,
                             methods.list.all,
                             weight = NULL,
                             stat.null = NULL,
                             null.max = 10^5,
                             Z.perm.all = NULL,
                             statistic = TRUE,
                             opt.method = "ILP_gurobi"){
  
  if(opt.method == "ILP_gurobi"){
    exact = TRUE
  } else {
    exact = FALSE
  }
  
  N = length(Y)
  p = N - k
  m = sum(Z)
  
  p_switch = p - 1
  
  
  if(!is.factor(block)){
    block = as.factor(block)
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
  
  coeflists_trt = comb_matrix_block(Z = Z, Y = Y,
                                    block = block, c = c,
                                    methods.list.all = methods.list.all)
  coeflists_con = comb_matrix_block(Z = 1 - Z, Y = -Y,
                                    block = block, c = c,
                                    methods.list.all = methods.list.all)
  
  #  switch_result = rep(NA, p_switch)
  pval = 1
  # p.trt = N - (N - m) - k
  # p.con = N - m - k
  
  for(iter in 0 : p_switch){
    if(iter > m || p_switch - iter > N- m) break    #condition for each pairs; (iter, N-k-1-iter)
    
    stat.min.trt = Gurobi_sol_com(Z = Z, block = block,
                                  weight = weight,
                                  coeflists = coeflists_trt,
                                  #p = p.trt,    ## should match with iterations?
                                  p = iter,   ## now k_1(and also k_0) differ for each iterations
                                  mu_sigma_list = ms_list,
                                  exact = exact)$obj
    pval.trt = mean(stat.null >= stat.min.trt)
    
    stat.min.control = Gurobi_sol_com(Z = 1 - Z, block = block,   ## should change this, since in optimization process Z is used for block treatment structure, which is switched in controlled case
                                      weight = weight,
                                      coeflists = coeflists_con,
                                      #p = p.con,
                                      p = p_switch - iter,
                                      mu_sigma_list = ms_list,
                                      exact = exact)$obj
    pval.control = mean(stat.null >= stat.min.control)
    
    stat.min.tmp = ifelse(pval.control < pval.trt, stat.min.trt, stat.min.control) ## take maximum
    pval.tmp = max(pval.trt, pval.control)
    
    stat.min = ifelse(pval < pval.tmp, stat.min, stat.min.tmp) ## take minimum
    pval = ifelse(pval < pval.tmp, pval, pval.tmp)
    print(iter)
  }
  pval = 2 * pval
  
  if(statistic == TRUE) {
    result = c(pval, stat.min)
    names(result) = c("p.value", "test.stat")
    return(result)
  } else return(pval)
}


#double check the function

N = 6
m = N / 2
#m = floor(150 * 3368 / 5768)

m.vec = floor(seq(2, m, len = 2))
alpha = 0.05
k = floor(N * 0.95)
c = 0

H = 2
s = 1

method.list.all.1 <- list()
method.list.all.2 <- list()

methods.list = list(method.lists.1,
                    method.lists.2)

for (j in 1 : H){
  for(i in 1 : s){
    methods.list[[j]][[i]] = list(name = "Polynomial",
                                  r = m.vec[j],
                                  std = TRUE,
                                  scale = TRUE)
  }
}



dumb.prop = 0.5 ## proportion of effected units
eff.len = floor(N * dumb.prop)    ### since k on tau_k is about the whole unit
eff.size = 1

eff.b = sample( rep( c(0, eff.size), c(N - eff.len, eff.len) ,
                     size = N), replace = FALSE)

Y0_tmp = rep(0, N)
Y1_tmp = (1 - eff.b) * Y0_tmp + eff.b

trt.b = sample(rep( c(0,1), c(N - m, m),
                    size = N), replace = FALSE)

Y = Y0_tmp * (1 - trt.b) + trt.b * Y1_tmp

Z_block = factor(rep(1, N))

switch_pval_comb(Z = trt.b,
                 Y = Y,
                 k = k,
                 c = 0,
                 block = Z_block,
                 methods.list.all = methods.list,
                 null.max = 10^2,
                 statistic = FALSE,
                 opt.method = "ILP_gurobi")



