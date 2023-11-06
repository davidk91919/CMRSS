
test_that("compare result with QIoT on SCRE", {
  
  # generating random data
  
  s = 5
  n = 6
  m = 0.5 * n
  N = s * n
  k = ceiling(0.95 * N)
  p = N - k
  c = 0
  
  simulation_Z <- as.data.frame(block.random(n = N, c(block = s))) 
  
  simulation_Z$block <- as.factor(simulation_Z$block)
  Z = block_ra(blocks = simulation_Z$block)
  Z_block <- simulation_Z$block
  
  Y0 = rep(NA, N)
  Y1 = rep(NA, N)
  for(i in 1 : s){
    Y0[Z_block == i] = rnorm(n)
    Y1[Z_block == i] = rnorm(n) +1
  }
  Y = Z * Y1 + (1-Z) * Y0
  
  # listing methods
  
  method.list.all.1 <- list()
  method.list.all.2 <- list()
  method.list.all.3 <- list()
  
  comb.methods.list.all <- list(method.list.all.1)
  
  methods.list.all <- list()
  
  H = 1 ## just one
  r.seq = ceiling(seq(2, n, len = H))
  
  
  for (j in 1 : H){
    for(i in 1 : s){
      comb.methods.list.all[[j]][[i]] = list(name = "Stephenson",
                                             s = r.seq[j],
                                             scale = FALSE)
    }
  } 
  for (j in 1 : H){
    for(i in 1 : s) {
      methods.list.all[[i]] <- list(name = "Stephenson",
                                    s = r.seq[j])
    }
  }
  
  expect_equal(
    pval_comb(Z = Z, Y = Y, k = k, c = c,
              block = Z_block,
              methods.list.all = comb.methods.list.all, 
              statistic = FALSE),
    
    pval_quantile_scre(Z = Z, Y = Y, block = Z_block, 
                       k = k, c = c, 
                       method.list.all = methods.list.all, opt.method = "ILP_gurobi")$lower,
    tolerance = 0.001)
  
})