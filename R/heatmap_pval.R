
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
                               Z.perm.all = Z.perm.all)
  }  
  colnames(res.mat) = round(c.vec, 2)
  rownames(res.mat) = k.vec
  
  pheatmap(res.mat,
           cluster_rows = FALSE,
           cluster_cols = FALSE)    ## omitted clustering. can add an option for this
}
