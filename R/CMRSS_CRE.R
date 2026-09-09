#' Sort treated units by outcome rank
#'
#' Returns the indices of treated units sorted by their outcome values
#' in increasing order.
#'
#' @param Y An n-dimensional observed outcome vector.
#' @param Z An n-dimensional binary treatment assignment vector.
#'
#' @return Integer vector of indices of treated units, sorted by increasing outcome.
#'
#' @examples
#' Z <- c(0, 0, 0, 1, 1, 1)
#' Y <- c(1, 2, 3, 30, 10, 20)
#' sort_treat(Y, Z)   # 5, 6, 4: the treated units ordered by outcome
#'
#' @export
sort_treat <- function(Y, Z){
  r <- rank(Y, ties.method = "first")
  ind.sort <- sort.int(r, index.return = TRUE)$ix
  ind.sort.treat <- ind.sort[Z[ind.sort] == 1]
  return(ind.sort.treat)
}



#' Generate complete randomization assignments
#'
#' Generates a matrix of permuted treatment assignments for completely
#' randomized experiments.
#'
#' @param n Total number of units.
#' @param m Number of treated units.
#' @param nperm Number of permutations. If Inf, generates all possible assignments.
#'
#' @return An n x nperm matrix where each column is a permuted treatment assignment.
#'
#' @export
assign_CRE <- function(n, m, nperm){
  if(is.finite(nperm)){
    # Vectorized: generate all sample indices at once
    indices <- replicate(nperm, sample.int(n, m, replace = FALSE))
    # indices is m x nperm matrix

    # Create Z.perm using matrix indexing
    Z.perm <- matrix(0, nrow = n, ncol = nperm)
    col_indices <- rep(1:nperm, each = m)
    row_indices <- as.vector(indices)
    Z.perm[cbind(row_indices, col_indices)] <- 1
  }

  if(is.infinite(nperm)){
    comb.all <- utils::combn(n, m)
    nperm <- ncol(comb.all)
    Z.perm <- matrix(0, nrow = n, ncol = nperm)
    # Vectorized using matrix indexing
    col_indices <- rep(1:nperm, each = m)
    row_indices <- as.vector(comb.all)
    Z.perm[cbind(row_indices, col_indices)] <- 1
  }

  return(Z.perm)
}



#' Calculate rank scores
#'
#' Computes rank scores for rank-based test statistics. Supports Wilcoxon,
#' Stephenson, and Polynomial rank scores.
#'
#' @param n Number of units.
#' @param method.list A list specifying the scoring method:
#'   \itemize{
#'     \item \code{name}: "Wilcoxon", "Stephenson", or "Polynomial"
#'     \item \code{s}: (for Stephenson) the parameter s
#'     \item \code{r}: (for Polynomial) the power parameter
#'     \item \code{std}: (for Polynomial) logical, use Puri(1965) normalization
#'     \item \code{scale}: logical, standardize scores to mean 0 and sd 1
#'   }
#'
#' @return A numeric vector of length n containing the rank scores.
#'
#' @details
#' Polynomial scores optionally use Puri(1965) normalization.
#' Stephenson scores use binomial coefficients.
#' Wilcoxon scores are simply the ranks 1 to n.
#'
#' @examples
#' rank_score(5, list(name = "Wilcoxon", scale = FALSE))
#' rank_score(5, list(name = "Stephenson", s = 3, scale = FALSE))
#' rank_score(5, list(name = "Polynomial", r = 2, std = TRUE, scale = FALSE))
#'
#' @export
rank_score <- function(n, method.list = list(name = "Polynomial", r, std = TRUE, scale = FALSE) ){
  if(method.list$name == "Polynomial"){
    r = method.list$r
    if(isTRUE(method.list$std)) {
      #      score = (c(1:n) / n)^(r-1)
      score = (c(1:n) / (n + 1) )^(r-1)  # changed
    } else {
      score = (c(1:n))^(r-1)
    }
    if(isTRUE(method.list$scale)){
      score = as.numeric(scale(score)) # normalized score; corresponds to 11p of the draft
    }
  }
  if(method.list$name == "Stephenson"){
    score = choose( c(1:n) - 1, method.list$s - 1)
    if(isTRUE(method.list$scale)){
      score = as.numeric(scale(score))
    }
  }
  if(method.list$name == "Wilcoxon"){
    score = c(1:n)
    if(isTRUE(method.list$scale)){
      score = as.numeric(scale(score))
    }
  }
  return(as.numeric(score))
}



#' Generate null distribution for rank sum statistic
#'
#' Computes the randomization null distribution of a single rank sum
#' statistic under complete randomization.
#'
#' @param n Total number of units.
#' @param m Number of treated units.
#' @param method.list A list specifying the scoring method.
#' @param score Optional pre-computed score vector.
#' @param nperm Number of permutations.
#' @param Z.perm Optional pre-computed permutation matrix.
#' @param chunk_size Integer chunk size used to generate permutations when
#'   \code{Z.perm} is \code{NULL}. Smaller values use less memory.
#'
#' @return A numeric vector of length nperm containing the null distribution.
#'
#' @examples
#' set.seed(1)
#' score <- rank_score(20, list(name = "Wilcoxon", scale = FALSE))
#' nd <- null_dist(n = 20, m = 12, score = score, nperm = 500)
#' quantile(nd, c(0.5, 0.95))
#'
#' @export
null_dist <- function(n, m, method.list = NULL, score = NULL,
                      nperm = 10^5, Z.perm = NULL, chunk_size = 1000){
  if(is.null(score)){
    score = rank_score( n, method.list )
  }
  score = as.numeric(score)

  if(!is.null(Z.perm)){
    return(as.vector(crossprod(score, Z.perm)))
  }

  if(!is.finite(nperm)){
    Z.perm = assign_CRE(n, m, nperm)
    return(as.vector(crossprod(score, Z.perm)))
  }

  if(m == 0){
    return(rep(0, nperm))
  }
  if(m == n){
    return(rep(sum(score), nperm))
  }

  stat.null = numeric(nperm)
  chunk_size = max(1L, as.integer(chunk_size))
  for(offset in seq.int(1L, nperm, by = chunk_size)){
    size = min(chunk_size, nperm - offset + 1L)
    indices = replicate(size, sample.int(n, m, replace = FALSE), simplify = "matrix")
    stat.null[offset:(offset + size - 1L)] = colSums(matrix(score[as.vector(indices)], nrow = m))
  }

  return(stat.null)
}

#' Generate null distributions for multiple rank sum statistics
#'
#' Computes the randomization null distribution for multiple rank sum
#' statistics simultaneously under complete randomization.
#'
#' @param n Total number of units.
#' @param m Number of treated units.
#' @param methods.list A list of method specifications.
#' @param nperm Number of permutations.
#' @param Z.perm Optional pre-computed permutation matrix.
#' @param chunk_size Integer chunk size used to generate permutations when
#'   \code{Z.perm} is \code{NULL}. Smaller values use less memory.
#'
#' @return An H x nperm matrix where H is the number of statistics.
#'
#' @examples
#' set.seed(1)
#' methods.list <- list(list(name = "Wilcoxon", scale = FALSE),
#'                      list(name = "Stephenson", s = 3, scale = FALSE))
#' dim(null_dist_multiple(n = 20, m = 12, methods.list = methods.list,
#'                        nperm = 500))
#'
#' @export
null_dist_multiple <- function(n, m, methods.list = NULL, nperm = 10^5, Z.perm = NULL, chunk_size = 1000){
  H = length(methods.list)

  if(!is.null(Z.perm)){
    nperm = ncol(Z.perm)
    stat.null.mult = matrix(NA, nrow = H, ncol = nperm)
    for(h in 1:H){
      stat.null.mult[h, ] = null_dist(n, m, method.list = methods.list[[h]], score = NULL, nperm = nperm, Z.perm = Z.perm)
    }
    return(stat.null.mult)
  }

  if(!is.finite(nperm)){
    Z.perm = assign_CRE(n, m, nperm)
    stat.null.mult = matrix(NA, nrow = H, ncol = ncol(Z.perm))
    for(h in 1:H){
      score = rank_score(n, methods.list[[h]])
      stat.null.mult[h, ] = as.vector(crossprod(score, Z.perm))
    }
    return(stat.null.mult)
  }

  if(m == 0){
    return(matrix(0, nrow = H, ncol = nperm))
  }
  if(m == n){
    result = matrix(NA, nrow = H, ncol = nperm)
    for(h in 1:H){
      score = as.numeric(rank_score(n, methods.list[[h]]))
      result[h, ] = rep(sum(score), nperm)
    }
    return(result)
  }

  scores = lapply(methods.list, function(method.list) rank_score(n, method.list))
  scores = lapply(scores, as.numeric)
  stat.null.mult = matrix(NA, nrow = H, ncol = nperm)

  chunk_size = max(1L, as.integer(chunk_size))
  for(offset in seq.int(1L, nperm, by = chunk_size)){
    size = min(chunk_size, nperm - offset + 1L)
    indices = replicate(size, sample.int(n, m, replace = FALSE), simplify = "matrix")
    for(h in 1:H){
      stat.null.mult[h, offset:(offset + size - 1L)] = colSums(matrix(scores[[h]][as.vector(indices)], nrow = m))
    }
  }

  return(stat.null.mult)
}



#' Compute minimum test statistic
#'
#' Calculates the minimum value of the rank sum statistic under the
#' null hypothesis, assuming at most n-k treated units can have effects
#' greater than c.
#'
#' @param Z An n-dimensional binary treatment assignment vector.
#' @param Y An n-dimensional observed outcome vector.
#' @param k Quantile index (between 1 and n).
#' @param c Threshold for the null hypothesis.
#' @param method.list A list specifying the scoring method.
#' @param score Optional pre-computed score vector.
#' @param ind.sort.treat Optional pre-computed sorted treatment indices.
#'
#' @return The minimum test statistic value.
#'
#' @references
#' Caughey, D., Dafoe, A., Li, X., & Miratrix, L. (2023).
#'
#' @examples
#' Z <- c(0, 0, 0, 1, 1, 1)
#' Y <- c(1, 2, 3, 10, 20, 30)
#' score <- rank_score(6, list(name = "Wilcoxon", scale = FALSE))
#' # k = 6 exempts no treated unit, so this is the treated rank sum: 4 + 5 + 6
#' min_stat(Z, Y, k = 6, c = 0, score = score)
#'
#' @export
min_stat <- function(Z, Y, k, c, method.list = NULL,
                     score = NULL,
                     ind.sort.treat = NULL){
  n = length(Y)
  m = sum(Z)

  # if(is.null(method.list)){   #added # I delete, we can have scale in the method.list
  #   score= scale(score)
  # }

  if(is.null(score)){
    score = rank_score(n, method.list)
  }

  # sort the treated units
  if(is.null(ind.sort.treat)){
    ind.sort.treat = sort_treat(Y, Z)
  }

  # get xi vector
  xi = rep(c, n)
  if(k < n){
    xi[ ind.sort.treat[ ( m + 1 - min(m, n-k) ):m ] ] = Inf
  }

  stat.min = sum( score[ rank( Y - Z * xi, ties.method = "first" )[Z==1] ] )

  return(stat.min)
}



#' Calculate p-value for single method
#'
#' Computes the p-value assuming at most n-k treated units have effects > c,
#' using a single rank sum statistic.
#'
#' @param Z An n-dimensional binary treatment assignment vector.
#' @param Y An n-dimensional observed outcome vector.
#' @param k Quantile index for the hypothesis.
#' @param c Threshold for the null hypothesis.
#' @param method.list A list specifying the rank score method.
#' @param score Optional pre-computed score vector.
#' @param stat.null Optional pre-computed null distribution.
#' @param nperm Number of permutations for null distribution.
#' @param Z.perm Optional permutation matrix.
#' @param ind.sort.treat Optional sorted indices of treated units.
#'
#' @return A numeric p-value.
#'
#' @keywords internal
pval_cre <- function(Z, Y, k, c,
                     method.list, score = NULL, stat.null = NULL,
                     nperm = 10^3, Z.perm = NULL, ind.sort.treat = NULL){
  n = length(Z)
  m = sum(Z)

  # get score if score is null #
  if(is.null(score)){
    score = rank_score( n, method.list )
  }

  # emp null dist #
  if(is.null(stat.null)){
    if(!is.null(Z.perm)){
      nperm = ncol(Z.perm)
    }
    stat.null = null_dist(n, m, score = score, nperm = nperm, Z.perm = Z.perm, method.list = method.list )
  }

  # min stat value under H_{k,c} #
  stat.min = min_stat(Z, Y, k, c, score = score, ind.sort.treat = ind.sort.treat, method.list = method.list )

  # p-value #
  pval = mean( stat.null >= stat.min )

  return(pval)
}



#' Null distribution for minimum p-values
#'
#' Computes the distribution of the minimum p-value, particularly Monte Carlo
#' samples from the combined null distribution F in Theorem 1 of the paper.
#'
#' @param n Total number of units.
#' @param m Number of treated units.
#' @param methods.list A list of method specifications.
#' @param Z.perm Optional permutation matrix.
#' @param nperm Number of permutations for null distribution.
#' @param stat.null.mult Optional pre-computed null distribution matrix.
#'
#' @return A numeric vector of minimum p-values under the null. Its length is
#'   the number of draws actually available: \code{ncol(Z.perm)} when a
#'   permutation matrix is supplied, \code{ncol(stat.null.mult)} when a null
#'   distribution is, and \code{nperm} otherwise.
#'
#' @examples
#' set.seed(1)
#' methods.list <- list(list(name = "Wilcoxon", scale = FALSE),
#'                      list(name = "Stephenson", s = 3, scale = FALSE))
#' cnd <- comb_null_dist_cre(n = 20, m = 12, methods.list = methods.list,
#'                           nperm = 500)
#' quantile(cnd, c(0.05, 0.10))
#'
#' @export
comb_null_dist_cre = function(n, m, methods.list, Z.perm = NULL, nperm = 10^4, stat.null.mult = NULL){

  H = length(methods.list)

  if(is.null(stat.null.mult)){
    stat.null.mult = null_dist_multiple(n, m, methods.list, nperm, Z.perm)
  }

  # How many draws there are is a property of the null distribution we now
  # hold, not of the nperm argument.  null_dist_multiple() already takes
  # nperm from ncol(Z.perm) when a permutation matrix is supplied, so reading
  # the width here keeps the two in step.  Before 0.2.10 a caller who passed
  # Z.perm together with a larger nperm got a tail.prob matrix of the wrong
  # width and R recycled the values into it, which made comb_p_val_cre return
  # p = 0 on data where the same permutations gave p = 0.8.
  nperm = ncol(stat.null.mult)

  tail.prob = matrix(NA, nrow = H, ncol = nperm)
  for(j in 1 : H){
    tail.prob[j, ] = 1 - (rank(stat.null.mult[j,], ties.method = "min") - 1) / nperm
  }

  result = apply(tail.prob, 2, min)
  return(result)
}



#' Calculate minimum p-value from multiple rank sum statistics
#'
#' Computes the minimum of p-values calibrated by the distribution of minimum
#' of tail probabilities.
#'
#' @param Z An n-dimensional binary treatment assignment vector.
#' @param Y An n-dimensional observed outcome vector.
#' @param k Quantile index for the hypothesis.
#' @param c Threshold for the null hypothesis.
#' @param methods.list A list of method specifications.
#' @param Z.perm Optional permutation matrix.
#' @param nperm Number of permutations for null distribution.
#' @param stat.null.mult Optional pre-computed null distribution matrix.
#'
#' @return A numeric minimum p-value.
#'
#' @keywords internal
min_p_multiple_rank_sum <- function(Z, Y, k, c, methods.list, Z.perm = NULL, nperm, stat.null.mult = NULL){
  
  n = length(Z)
  m = sum(Z)
  
  if(is.null(stat.null.mult)){
    stat.null.mult = null_dist_multiple(n, m, methods.list, nperm, Z.perm)
  }

  H = length(methods.list)
  pval.vec = rep(NA, H)
  for (i in 1 : H){
    method.list = methods.list[[i]]
    pval.vec[i] = pval_cre(Z = Z, Y = Y, k = k, c = c, Z.perm = Z.perm, method.list = method.list, stat.null = stat.null.mult[i, ])
  }

  min.pval = min(pval.vec)
  return(min.pval)
}


#' Combined p-value for quantile treatment effects in CRE
#'
#' Calculate a valid p-value, based on multiple rank sum statistics,
#' for testing the null hypothesis about quantiles of individual treatment
#' effects in completely randomized experiments (CRE).
#'
#' @description
#' Tests the null hypothesis \eqn{H_0: \tau_{(k)} \leq c}, where \eqn{\tau_{(k)}}
#' denotes the individual treatment effect at rank k. The test combines multiple
#' rank sum statistics to improve power.
#'
#' @param Z An n-dimensional binary treatment assignment vector (1 = treated, 0 = control).
#' @param Y An n-dimensional observed outcome vector.
#' @param k An integer between 1 and n specifying which quantile of
#'   individual effect is of interest.
#' @param c A numeric value specifying the threshold for the null hypothesis.
#' @param methods.list A list of method specifications for the rank sum statistics.
#'   Each element should be a list with:
#'   \itemize{
#'     \item \code{name}: "Wilcoxon", "Stephenson", or "Polynomial"
#'     \item \code{s}: (for Stephenson) parameter s
#'     \item \code{r}: (for Polynomial) power parameter
#'     \item \code{std}: (for Polynomial) logical, use Puri normalization
#'     \item \code{scale}: logical, standardize scores
#'   }
#' @param Z.perm An n x nperm matrix of permuted treatment assignments for
#'   approximating the null distribution. If NULL, generated automatically.
#' @param nperm Number of permutations for approximating the null distribution.
#' @param stat.null.mult A matrix whose empirical distribution approximates
#'   the randomization distribution of multiple rank statistics. If NULL,
#'   computed from Z.perm.
#'
#' @return A numeric p-value for testing the specified null hypothesis.
#'
#' @examples
#' # Simple example with Wilcoxon and Stephenson statistics
#' set.seed(123)
#' n <- 30
#' Z <- sample(c(rep(1, 15), rep(0, 15)))
#' Y <- rnorm(n) + 0.5 * Z  # Treatment effect of 0.5
#'
#' # Define methods: Wilcoxon and Stephenson with s=3
#' methods.list <- list(
#'   list(name = "Wilcoxon", scale = FALSE),
#'   list(name = "Stephenson", s = 3, scale = FALSE)
#' )
#'
#' # Test if the 80th percentile effect is <= 0
#' k <- floor(0.8 * n)
#' pval <- comb_p_val_cre(Z, Y, k, c = 0, methods.list, nperm = 1000)
#'
#' @seealso \code{\link{pval_comb_block}} for stratified experiments
#' @export
comb_p_val_cre = function(Z, Y, k, c, methods.list,
                          Z.perm = NULL, nperm,
                          stat.null.mult = NULL
                          ){
  H = length(methods.list)
  pval.vec = rep(NA, H)
  n = length(Z)
  m = sum(Z)

  if(is.null(stat.null.mult)){
    stat.null.mult = null_dist_multiple(n, m, methods.list, nperm, Z.perm)
  }

  null.dist.comb = comb_null_dist_cre(n, m, methods.list, Z.perm, nperm, stat.null.mult)

  min.pval = min_p_multiple_rank_sum(Z, Y, k, c, methods.list, Z.perm, nperm, stat.null.mult)
  result = mean(null.dist.comb <= min.pval) # XL: I changed
  return(result)
}


#' Simultaneous inference for multiple quantiles on CRE (treated units)
#'
#' Helper function for computing lower limits of prediction intervals for
#' quantiles across the treated units using combined p-value.
#'
#' @param Z An n-dimensional binary treatment assignment vector.
#' @param Y An n-dimensional observed outcome vector.
#' @param methods.list A list of method specifications.
#' @param nperm Number of permutations for null distribution.
#' @param k.vec Vector of quantile indices to compute intervals for.
#' @param Z.perm Optional permutation matrix.
#' @param alpha Significance level.
#' @param tol Tolerance for root-finding.
#' @param ind.sort.treat Optional sorted indices of treated units.
#'
#' @return A numeric vector of lower confidence limits.
#'
#' @keywords internal
com_conf_quant_larger_trt <- function( Z, Y, methods.list = NULL,
                                       nperm = 10^4,
                                       k.vec = NULL,
                                       Z.perm = NULL,
                                       alpha = 0.05,
                                       tol = 10^(-3),
                                       ind.sort.treat = NULL ){
  n = length(Z)
  m = sum(Z)

  # emp null dist #
  # XL: I add the following
  if(!is.null(Z.perm)){
    nperm = ncol(Z.perm)
  }

  stat.null.mult = null_dist_multiple(n, m, methods.list, nperm, Z.perm)
  null.dist.comb = comb_null_dist_cre(n, m, methods.list, Z.perm, nperm, stat.null.mult)

  # find threshold such that valid p-value <= alpha <===> min_pval < something <===>  -min_pval > -something #
  thres.minus.min.pval = sort(-1*null.dist.comb, decreasing = TRUE)[ floor(nperm * alpha) + 1 ]
  thres = alpha

  # range of c #
  Y1.max = max(Y[Z==1])
  Y1.min = min(Y[Z==1])
  Y0.max = max(Y[Z==0])
  Y0.min = min(Y[Z==0])

  c.max = Y1.max - Y0.min + tol
  c.min = Y1.min - Y0.max - tol


  # sort the treated units
  if(is.null(ind.sort.treat)){
    ind.sort.treat = sort_treat(Y, Z)
  }

  if( is.null(k.vec)) {
    # conf interval for all quantiles tau_{(k)} #
    c.limit = rep(NA, n)

    for(k in n:(n-m)){
      # define the target fun #
      # f > 0 <==> p-value <= alpha
      # f decreases in c, p value increases in c
      # we want to find minimum c which is close as possible to f <= 0; for 1-alpha CI
      if(k < n){
        c.max = c.limit[k+1]
      }
      ## XL: an alternative way that using minus min pval as the test statistic
      f <- function(c){
        stat.min = -1 * min_p_multiple_rank_sum(Z, Y, k, c, methods.list, Z.perm, nperm, stat.null.mult)
        return(stat.min - thres.minus.min.pval)
      }

      # check whether f(-Inf) = f(c.min) > 0 #
      if( f(c.min) <= 0 ){
        c.sol = -Inf
      }
      if( f(c.max) > 0){
        c.sol = c.max
      }
      if( f(c.min) > 0 & f(c.max) <= 0 ){
        c.sol = uniroot(f, interval = c(c.min, c.max), extendInt = "downX", tol = tol)$root
        # find the min c st p-value > alpha <==> f <= 0 #
        c.sol = round(c.sol, digits = -log10(tol))
        if( f(c.sol) <= 0 ){
          while( f(c.sol) <= 0 ){
            c.sol = c.sol - tol
          }
          c.sol = c.sol + tol
        }else{
          while(f(c.sol) > 0 ){
            c.sol = c.sol + tol
          }
        }
      }

      c.limit[k] = c.sol
    }

    if( n-m > 1){
      c.limit[1:(n-m-1)] = c.limit[n-m]
    }
  }else{

    k.vec.sort = sort(k.vec, decreasing = FALSE)
    j.max = length(k.vec.sort)
    j.min = max( sum(k.vec <= (n-m)), 1)

    c.limit = rep(NA, j.max)

    for(j in j.max:j.min){
      k = k.vec.sort[j]
      if(j < j.max){
        c.max = c.limit[j+1]
      }
      ## XL: an alternative way that using minus min pval as the test statistic
      f <- function(c){
        stat.min = -1 * min_p_multiple_rank_sum(Z, Y, k, c, methods.list, Z.perm, nperm, stat.null.mult)
        return(stat.min - thres.minus.min.pval)
      }

      # check whether f(-Inf) = f(c.min) > 0 #
      if( f(c.min) <= 0 ){
        c.sol = -Inf
      }
      if( f(c.max) > 0){
        c.sol = c.max
      }
      if( f(c.min) > 0 & f(c.max) <= 0 ){
        c.sol = uniroot(f, interval = c(c.min, c.max), extendInt = "downX", tol = tol)$root
        # find the min c st p-value > alpha <==> f <= 0 #
        c.sol = round(c.sol, digits = -log10(tol))
        if( f(c.sol) <= 0 ){
          while( f(c.sol) <= 0 ){
            c.sol = c.sol - tol
          }
          c.sol = c.sol + tol
        }else{
          while(f(c.sol) > 0 ){
            c.sol = c.sol + tol
          }
        }
      }

      c.limit[j] = c.sol
    }
    if(j.min > 1){
      c.limit[1:(j.min-1)] = c.limit[j.min]
    }
  }

  return( c.limit )
}


#' Simultaneous lower bounds using multiple rank sum statistics in CRE
#'
#' Computes simultaneous confidence/prediction intervals for quantiles of
#' individual treatment effects in completely randomized experiments (CRE)
#' by combining multiple rank sum statistics.
#'
#' @param Z An \eqn{n} dimensional treatment assignment vector (1 = treated, 0 = control).
#' @param Y An \eqn{n} dimensional observed outcome vector.
#' @param methods.list A list of lists specifying the choice of multiple
#'   rank sum test statistics. Each element should be a list with:
#'   \itemize{
#'     \item \code{name}: "Wilcoxon", "Stephenson", or "Polynomial"
#'     \item \code{s}: (for Stephenson) parameter s controlling sensitivity to upper ranks
#'     \item \code{r}: (for Polynomial) power parameter
#'     \item \code{std}: (for Polynomial) logical, use Puri(1965) normalization
#'     \item \code{scale}: logical, standardize scores to mean 0 and sd 1
#'   }
#' @param nperm A positive integer representing the number of
#'   permutations for approximating the randomization distribution of
#'   the rank sum statistic.
#' @param set Set of quantiles of interest:
#'   \itemize{
#'     \item "treat": Prediction intervals for effect quantiles among treated units
#'     \item "control": Prediction intervals for effect quantiles among control units
#'     \item "all": Confidence intervals for all effect quantiles
#'   }
#' @param Z.perm A \eqn{n \times nperm} matrix that specifies the
#'   permuted assignments for approximating the null distribution of
#'   the test statistic. If NULL, generated automatically.
#' @param alpha A numeric value where 1-alpha indicates the
#'   confidence level.
#' @param tol A numeric value specifying the precision of the
#'   obtained confidence intervals. For example, if tol = 10^(-3),
#'   then the confidence limits are precise up to 3 digits.
#'
#' @return A vector specifying lower limits of prediction (confidence) intervals for
#'   quantiles k = 1 ~ m (for "treat"), n - m + 1 ~ n (for "control"), or 1 ~ n (for "all").
#'
#' @details
#' This function implements the combined rank sum test approach for inference
#' about quantiles of individual treatment effects. By combining multiple rank
#' statistics (e.g., Stephenson statistics with different s values), the method
#' can achieve better power across a range of effect distributions.
#'
#' When \code{set = "all"}, the function combines inference from both treated and
#' control units using the approach described in Chen and Li (2024), with a
#' Bonferroni-style adjustment (alpha/2 for each direction).
#'
#' @examples
#' \dontrun{
#' # Load the electric teachers dataset
#' data(electric_teachers)
#'
#' # Set up treatment and outcome (treating as CRE, ignoring sites)
#' Z <- electric_teachers$TxAny
#' Y <- electric_teachers$gain
#'
#' # Define multiple Stephenson statistics with different s values
#' # Larger s focuses more on upper ranks (larger treatment effects)
#' s.vec <- c(2, 6, 10, 30)
#' methods.list <- lapply(s.vec, function(s) {
#'   list(name = "Stephenson", s = s, std = TRUE, scale = TRUE)
#' })
#'
#' # Prediction intervals for treated units (90% confidence)
#' ci.treat <- com_conf_quant_larger_cre(Z, Y,
#'                                       methods.list = methods.list,
#'                                       nperm = 10000,
#'                                       set = "treat",
#'                                       alpha = 0.05)
#'
#' # Prediction intervals for control units
#' ci.control <- com_conf_quant_larger_cre(Z, Y,
#'                                         methods.list = methods.list,
#'                                         nperm = 10000,
#'                                         set = "control",
#'                                         alpha = 0.05)
#'
#' # Confidence intervals for all effect quantiles
#' ci.all <- com_conf_quant_larger_cre(Z, Y,
#'                                     methods.list = methods.list,
#'                                     nperm = 10000,
#'                                     set = "all",
#'                                     alpha = 0.10)
#' }
#'
#' @seealso \code{\link{com_block_conf_quant_larger}} for stratified experiments,
#'   \code{\link{comb_p_val_cre}} for p-value computation
#' @export
com_conf_quant_larger_cre <- function( Z, Y, methods.list,
                                       nperm = 10^4,
                                       set = "treat",
                                       Z.perm = NULL,
                                       alpha = 0.05,
                                       tol = 10^(-3)
){

  H = length(methods.list)
  n = length(Z)

  # prediction intervals for effect quantiles among treated units
  if(set == "treat"){
    ci.treat = com_conf_quant_larger_trt(Z, Y, methods.list,
                                         nperm,
                                         k.vec = NULL,
                                         Z.perm, alpha, tol, ind.sort.treat = NULL
    )
    ci.treat = ci.treat[(n - sum(Z) + 1):n]
    return(ci.treat)
  }

  # prediction intervals for effect quantiles among control units
  if(set == "control"){
    # stat.null.control = NULL
    if(is.null(Z.perm)){
      Z.perm.control = NULL
    }else{
      Z.perm.control = 1 - Z.perm
    }

    Z = 1 - Z
    Y = -Y
    ci.control = com_conf_quant_larger_trt(Z, Y,
                                           methods.list = methods.list,
                                           nperm = nperm,
                                           k.vec = NULL,
                                           Z.perm = Z.perm.control,
                                           alpha = alpha, tol = tol, ind.sort.treat = NULL)
    ci.control = ci.control[(n - sum(Z) + 1):n]
    return(ci.control)
  }

  # set = "all": with Chen and Li's adjustment
  if(set == "all"){
    # trted ones
    ci.treat = com_conf_quant_larger_trt(Z, Y, methods.list,
                                         nperm, k.vec = NULL,
                                         Z.perm, alpha = alpha/2, tol, ind.sort.treat=NULL)
    ci.treat = ci.treat[(n - sum(Z) + 1):n]

    # controlled ones
    # stat.null.control = NULL

    if(is.null(Z.perm)){
      Z.perm.control = NULL
    }else{
      Z.perm.control = 1 - Z.perm
    }
    Z = 1 - Z
    Y = -Y
    ci.control = com_conf_quant_larger_trt(Z, Y,
                                           methods.list = methods.list,
                                           nperm = nperm,
                                           k.vec = NULL,
                                           Z.perm = Z.perm.control,
                                           alpha = alpha/2, tol = tol, ind.sort.treat=NULL)
    ci.control = ci.control[(n - sum(Z) + 1):n]

    ci.all = sort( c(ci.treat, ci.control))
    return(ci.all)
  }
}
