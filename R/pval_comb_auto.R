#' pval_comb_auto
#'
#' Function to calculate desired p-value combining 10 methods without specification
#'
#' @export


pval_comb_auto = function(Z, Y, k, c,
                          block,
                          null.max = 10^3,
                          STRA_CODE,
                          statistic = FALSE) {

  if(any(is.na(Y))){stop("Outcome Should Not Include NA Values")}

  Y.obs = Y
  Z.obs = Z
  block.obs = factor(block)

  m = max(table(STRA_CODE, block, exclude=c())[2,])

  H = 10 # will combine 10 many rank sum statistics, using polynomial rank score
  s = length(levels(block.obs)))
  m.seq = ceiling(seq(2, m, len = H))

  # listing methods

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

  ## implement method

  return(pval_comb(Z = Z.obs,
                   Y = Y.obs,
                   k = k,
                   c = 0,
                   block = block.obs,
                   methods.list.all = methods.list,
                   null.max = 10^3,
                   statistic = FALSE))

}


