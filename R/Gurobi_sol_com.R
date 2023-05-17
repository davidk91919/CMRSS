#' Function implementing optimization for multiple stratified rank sum statistics, using Gurobi software.
#'
#' Obtain the value of test statistic by optimization, for \eqn{\underset{h}{max} \frac{t^{(h)}(Z, Y) - \mu_h }{\sigma_h}}.
#'
#' @param Z n-dimension treatment assignment vector, for every units.
#' @param block n-dimension vector specifies stratum status of every units.
#' @param weight Vector which has same length with total number of stratum. Giving weight on each stratum before optimization. Basic option is number of treated units on each stratum.
#' @param coeflists All possible test statistics value, which will usually be a output from \emph{comb_matrix_block} function.
#' @param p N-k, where k is given in null hypothesis as \eqn{H_{N,k,c} = \tau_{k} \leq c}.
#' @param mu_sigma_list List of population mean / variance of every stratified rank sum statistics. This will usually be a output from \emph{mu_sigma_list} function.
#' @param exact Option for type of optimization problem. If exact = TRUE, will try ILP(integer linear programming). If exact = FALSE, will relax the integer constraints and try linear programming. 
#'
#' @export

Gurobi_sol_com <- function(Z, block, weight = NULL, coeflists, p, mu_sigma_list, exact = TRUE){

  ## n: number of all units + B
  ## B: number of strata
  ## H: number of stat to combine
  ## Q is organized by x's(n*H), eta's(H), theta
  ## weight is a H by B matrix, where each hth row is a weight vector for hth statistic.

  model = list()
  H = length(coeflists)
  B = length(coeflists[[1]])
  nb = vector(length = B)

  coeflist = coeflists[[1]]

  for (i in 1 : B){
    nb[i] = ncol(coeflist[[i]])
  }

  n = sum(nb)
  indx = c(0, cumsum(rep(nb, length(coeflists))))
  Q = rep(0, n*H + H + 1)
  Q[n*H + H + 1] = 1

  block.levels = levels(block)
  ## calculate numbers of treated units(m_s) in each block
  mb = rep(NA, B)
  for(i in 1:B){
    mb[i] = sum(Z[block == block.levels[i]])
  }

  if(is.null(weight)){
    weight = matrix(nrow = H, ncol = B)
    for(i in 1 : H){
      weight[i,] = mb
    }
  }

  ## sum of x in each strata is 1
  Ai = rep(1:(H*B), rep(nb, H))
  Aj = 1:(n*H)
  x = rep(1, n*H)

  ## cost is N-k for all choices of stat and profit for h's stat is eta'h
  Ai = c(Ai, (H*(B+1)+1):(H*(B+2)))
  Aj = c(Aj, (n*H + 1):(n*H+H))
  x = c(x, rep(-1, H))
  for (h in 1 : H){
    coeflist = coeflists[[h]]
    for (i in 1 : B) {
      Ai = c(Ai, rep(c(H*B + h, H*(B+1)+h), each = indx[(h-1)*B+i+1]-indx[(h-1)*B+i]))
      Aj = c(Aj, rep((indx[(h-1)*B+i]+1): indx[(h-1)*B+i+1], 2))
      x = c(x, coeflist[[i]][1,], weight[h,i] * (1 / mb[i]) * coeflist[[i]][2,] )
    }
  }

  ## for final constraints
  mu = mu_sigma_list$mean
  sig_rev = 1/mu_sigma_list$sigma
  Ai = c(Ai, rep((H*(B+2)+1):(H*(B+3)),2))
  Aj = c(Aj, rep(n*H + H + 1, H), (n*H+1):(n*H+H))
  x = c(x, rep(-1, H), sig_rev)

  A = sparseMatrix(Ai, Aj, x = x)


  model$A = A
  model$obj = Q
  model$modelsense = "min"
  model$rhs = c(rep(1, H*B), rep(p, H), rep(0, H), mu*sig_rev)
  model$sense = c(rep("=", H*B), rep("<=", H), rep("=", H), rep("<=", H))
  model$lb = c(rep(0, n*H), rep(-Inf, H + 1))
  if(exact){
    model$vtype = c(rep("B", n*H), rep("C", (H + 1) ))
  }

  params <- list(OutputFlag = 0)
  result = gurobi::gurobi(model,params)

  return(list(sol = result$x, obj = result$objval ))

}
