
test_that("pvalue in one stratum", {
  Z = sample(c(0,1), sample(100, 1), replace = TRUE)
  n = length(Z)
  m = sum(Z)
  
  Y = sample(1:100, n, replace = TRUE)
  k = sample(n, 1)
  c = sample(c(-1,0,1), 1)
  
  r = sample(m, 1)
  block = factor(rep(1, n))
  
  method.list.old.wil = list(name = "Wilcoxon")
  method.list.old.ste = list(name = "Stephenson",
                             s = r)
  
  method.list.new.wil = list(name = "Wilcoxon", scale = FALSE)
  method.list.new.ste = list(name = "Stephenson", scale = FALSE,
                             s = r)
  
  method.lists.wil = list(list(method.list.new.wil))
  method.lists.ste = list(list(method.list.new.ste))
  
  expect_equal(pval_quantile(Z = Z, Y = Y, k = k, c = c, method.list = method.list.old.wil), 
               pval_comb(Z = Z, Y = Y, k = k, c = c,
                         block = block, methods.list.all = method.lists.wil, statistic = FALSE), tolerance = 0.001)
  
  expect_equal(pval_quantile(Z = Z, Y = Y, k = k, c = c, method.list = method.list.old.ste), 
               pval_comb(Z = Z, Y = Y, k = k, c = c,
                         block = block, methods.list.all = method.lists.ste, statistic = FALSE), tolerance = 0.001)
})