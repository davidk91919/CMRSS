#' heatmap_pval
#'
#' A function to visualize multiple p-values corresponding to multiple choices of upper bounds and quantiles.
#'
#' @export





### Function for plotting heatmap, corresponding to multiple k(quantile of individual treatment effect) and c(upper bound on the null)
heatmap_pval = function(Z, Y,
                        k.vec = NULL,
                        c.vec = NULL,
                        block,
                        methods.list.all,
                        auto = TRUE,
                        STRA_CODE,
                        quantile = FALSE,
                        weight = NULL,
                        stat.null = NULL,
                        null.max = 10^5,
                        Z.perm.all = NULL){


  if(!is.factor(block)){
    block = as.factor(block)
  }

  if(is.null(k.vec)){
    k.vec = floor(length(Y) * seq(0, 1, len = 5))
  }
  if(is.null(c.vec)){
    c.vec = -2:2
  }

  l = length(k.vec)
  v = length(c.vec)

  res.mat = data.frame(0, nrow = l, ncol = v)

  if(auto == TRUE){
    for(i in 1 : l){
      for(k in 1 : v){
    res.mat[i,k] = pval_comb_auto(Z = Z,
                                  Y = Y,
                                  k = k.vec[i],
                                  c = c.vec[k],
                                  block = block,
                                  STRA_CODE = STRA_CODE,
                                  statistic = FALSE)
      }
    }
  } else {
  for(i in 1 : l){
    for(k in 1 : v){
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
  }
  }
  res.mat = t(res.mat)
  rownames(res.mat) = round(c.vec, 2)
  if(quantile == FALSE){
  colnames(res.mat) = k.vec
  } else {
  colnames(res.mat) = round(k.vec / length(Y), 2)
  }


  pheatmap::pheatmap(res.mat,
           cluster_rows = FALSE,
           cluster_cols = FALSE)
}
