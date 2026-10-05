#' Summarize block structure
#'
#' Computes summary statistics for a stratified randomized experiment,
#' including the number of blocks, units per block, and treated units per block.
#'
#' @param Z An n-dimensional binary treatment assignment vector.
#' @param block An n-dimensional vector or factor specifying block membership.
#'
#' @return A list with components:
#'   \item{block}{The block factor}
#'   \item{B}{Number of blocks}
#'   \item{nb}{Vector of block sizes}
#'   \item{mb}{Vector of treated units per block}
#'   \item{mb_ctrl}{Vector of control units per block}
#'   \item{units.block}{List of unit indices per block}
#'   \item{block.levels}{Levels of the block factor}
#'
#' @export
summary_block <- function(Z, block){
  # number of blocks
  if(!is.factor(block)){
    block = as.factor(block)
  }
  block.levels = levels(block)
  B = length(block.levels)

  # number of treated and control within each block
  units.block = split(seq_along(block), block)
  nb = unname(vapply(units.block, length, integer(1)))
  mb = unname(vapply(units.block, function(idx) sum(Z[idx]), numeric(1)))
  units.block = unname(units.block)
  mb_ctrl = nb - mb

  result = list(block=block, B = B, nb = nb, mb = mb, mb_ctrl = mb_ctrl, units.block = units.block, block.levels = block.levels )

  return(result)
}

#' Generate block-randomized treatment assignments
#'
#' Generates permuted treatment assignments that respect the block structure
#' of a stratified randomized experiment.
#'
#' @param block.sum A block summary object from \code{summary_block}.
#' @param null.max Number of permutations to generate.
#'
#' @return An n x null.max matrix of permuted treatment assignments.
#'
#' @export
assign_block <- function(block.sum, null.max = 10^4){
  block = block.sum$block
  B = block.sum$B
  nb = block.sum$nb
  mb = block.sum$mb
  mb_ctrl = block.sum$mb_ctrl
  units.block = block.sum$units.block
  block.levels = block.sum$block.levels

  n = sum(nb)

  Z.perm = matrix(0, nrow = n, ncol = null.max)

  # Vectorized: for each block, generate all null.max permutations at once
  for(i in 1:B){
    block_vec = c(rep(1, mb[i]), rep(0, mb_ctrl[i]))
    # Generate null.max permutations for this block (nb[i] x null.max matrix)
    block_perms = replicate(null.max, sample(block_vec))
    # Assign to the appropriate rows
    Z.perm[units.block[[i]], ] = block_perms
  }

  return(Z.perm)
}

# Internal helper: generate a chunk of block-randomized assignments without
# materializing the full n x null.max permutation matrix.
.generate_Z_chunk_block <- function(units.block, mb, n, size, combo_cache = NULL, max_ncomb = 1e5){
  B = length(units.block)
  Z.chunk = matrix(0L, nrow = n, ncol = size)

  if(is.null(combo_cache)){
    combo_cache = new.env(parent = emptyenv())
  }

  for(i in 1:B){
    units = units.block[[i]]
    nbi = length(units)
    mi = as.integer(mb[i])

    if(mi <= 0L){
      next
    }
    if(mi >= nbi){
      Z.chunk[units, ] = 1L
      next
    }

    if(mi == 1L){
      treated = sample.int(nbi, size, replace = TRUE)
      Z.chunk[cbind(units[treated], seq_len(size))] = 1L
      next
    }
    if(mi == (nbi - 1L)){
      Z.chunk[units, ] = 1L
      ctrl = sample.int(nbi, size, replace = TRUE)
      Z.chunk[cbind(units[ctrl], seq_len(size))] = 0L
      next
    }

    key = paste0(nbi, ":", mi)
    comb = NULL
    if(exists(key, envir = combo_cache, inherits = FALSE)){
      comb = get(key, envir = combo_cache, inherits = FALSE)
    } else {
      log_ncomb = lchoose(nbi, mi)
      if(is.finite(log_ncomb) && log_ncomb <= log(max_ncomb)){
        comb = utils::combn(nbi, mi)
        assign(key, comb, envir = combo_cache)
      }
    }

    if(!is.null(comb)){
      ncomb = ncol(comb)
      picks = sample.int(ncomb, size, replace = TRUE)
      idx = comb[, picks, drop = FALSE]
    } else {
      idx = replicate(size, sample.int(nbi, mi, replace = FALSE), simplify = "matrix")
    }

    row_idx = units[as.vector(idx)]
    col_idx = rep(seq_len(size), each = mi)
    Z.chunk[cbind(row_idx, col_idx)] = 1L
  }

  return(Z.chunk)
}

#' Compute block weights
#'
#' Computes weights for combining block-specific statistics according
#' to a specified weighting scheme.
#'
#' @param block.sum A block summary object from \code{summary_block}.
#' @param weight.name Weighting scheme: "asymp.opt" for asymptotically optimal
#'   weights, or "dis.free" for distribution-free weights.
#'
#' @return A B-dimensional vector of block weights.
#'
#' @keywords internal
weight_scheme <- function(block.sum, weight.name = "asymp.opt"){
  B = block.sum$B
  nb = block.sum$nb
  mb = block.sum$mb
  mb_ctrl = block.sum$mb_ctrl

  weight = rep(NA, B)

  if(weight.name == "asymp.opt"){
    weight = mb
  }
  if(weight.name == "dis.free"){
    weight = (nb + 1) / mb_ctrl
  }

  return(weight)
}

#' Compute scores for all blocks
#'
#' Computes rank scores for each block according to block-specific methods.
#'
#' @param nb Vector of block sizes.
#' @param method.list.all A list of method specifications, one per block.
#'
#' @return A list of score vectors, one per block.
#'
#' @keywords internal
score_all_blocks <- function(nb, method.list.all){
  B = length(method.list.all)
  score.list.all = list()
  for(i in 1 : B){
    method.list = method.list.all[[i]]
    score = rank_score(nb[i], method.list)
    score.list.all[[i]] = score
  }
  return(score.list.all)
}


#' Compute weighted stratified rank sum statistic
#'
#' Calculates a single weighted stratified rank sum statistic combining
#' block-specific statistics.
#'
#' @param Z An n-dimensional binary treatment assignment vector.
#' @param Y An n-dimensional observed outcome vector.
#' @param block An n-dimensional vector specifying block membership.
#' @param method.list.all A list of method specifications, one per block.
#' @param score.list.all Optional pre-computed scores.
#' @param weight A B-dimensional vector of block weights.
#' @param block.sum Optional pre-computed block summary.
#'
#' @return A scalar test statistic value.
#'
#' @keywords internal
single_weight_strat_rank_sum_stat <- function(Z, Y, block, method.list.all, score.list.all = NULL, weight, block.sum = NULL){

  if(is.null(block.sum)){
    block.sum = summary_block(Z, block)
  }
  block = block.sum$block
  B = block.sum$B
  nb = block.sum$nb
  mb = block.sum$mb
  mb_ctrl = block.sum$mb_ctrl
  units.block = block.sum$units.block
  block.levels = block.sum$block.levels

  if(is.null(score.list.all)){
    score.list.all = score_all_blocks(nb, method.list.all)
  }

  result = 0
  for(i in 1 : B){
    score = score.list.all[[i]]


    zs = Z[units.block[[i]]]
    ys = Y[units.block[[i]]]

    stat.block = sum( score[rank(ys, ties.method = "first")[zs == 1] ] )
    result = result + weight[i] * 1 / mb[i] * stat.block
  }
  return(result)
}

#' Compute mean and standard deviation of stratified statistic
#'
#' Calculates the mean and standard deviation of a weighted stratified
#' rank sum statistic under the null hypothesis.
#'
#' @param Z An n-dimensional binary treatment assignment vector.
#' @param block An n-dimensional vector specifying block membership.
#' @param method.list.all A list of method specifications, one per block.
#' @param score.list.all Optional pre-computed scores.
#' @param weight A B-dimensional vector of block weights.
#' @param block.sum Optional pre-computed block summary.
#'
#' @return A list with components:
#'   \item{mu}{Overall weighted mean}
#'   \item{sigma}{Overall weighted standard deviation}
#'   \item{mu_weight}{Block-wise weighted means (mean_block * weight)}
#'   \item{mu_block}{Block-wise unweighted means}
#'   \item{sd_block}{Block-wise weighted standard deviations}
#'
#' @keywords internal
mu_sigma_single_weight_strat_rank_sum_stat <- function(Z, block,
                                                       method.list.all,
                                                       score.list.all = NULL,
                                                       weight,
                                                       block.sum = NULL){

  if(is.null(block.sum)){
    block.sum = summary_block(Z, block)
  }
  block = block.sum$block
  B = block.sum$B
  nb = block.sum$nb
  mb = block.sum$mb
  mb_ctrl = block.sum$mb_ctrl
  units.block = block.sum$units.block
  block.levels = block.sum$block.levels

  if(is.null(score.list.all)){
    score.list.all = score_all_blocks(nb, method.list.all)
  }


  mean_block = rep(NA, B)
  var_block = rep(NA, B)

  for(b in 1 : B){
    method.list = method.list.all[[b]]
    ns = nb[b]
    ms = mb[b]
    cs = mb_ctrl[b]

    score_b = score.list.all[[b]]
    mean_block[b] = mean(score_b)
    var_block[b] = cs/(ms*ns) * var(score_b)
  }

  mu_block = mean_block
  mu_weight = mean_block * weight
  single_mu = sum(mean_block * weight)

  single_sigma = sqrt(sum(weight^2 * var_block))
  sigma_block = sqrt( weight^2 * var_block)

  return( list(mu = single_mu,
               sigma = single_sigma,
               mu_weight = mu_weight,
               mu_block = mu_block,
               sd_block = sigma_block) )
}


#' Generate null distribution of combined stratified statistics
#'
#' Generates the randomization null distribution of the maximum among
#' weighted and normalized stratified rank sum statistics.
#'
#' @param Z An n-dimensional binary treatment assignment vector.
#' @param block An n-dimensional vector specifying block membership.
#' @param methods.list.all A list of lists of method specifications.
#' @param scores.list.all Optional pre-computed scores.
#' @param null.max Number of permutations.
#' @param weight A B-dimensional vector of block weights.
#' @param block.sum Optional pre-computed block summary.
#' @param Z.perm Optional pre-computed permutation matrix.
#' @param chunk_size Optional number of permutations per chunk when
#'   \code{Z.perm} is \code{NULL}. If \code{NULL}, chosen adaptively based on \code{n}.
#'
#' @return A numeric vector of the null distribution (length null.max).
#'
#' @keywords internal
com_null_dist_block <- function(Z, block,
                                methods.list.all,
                                scores.list.all = NULL,
                                null.max = 10^4,
                                weight,
                                block.sum = NULL,
                                Z.perm = NULL,
                                chunk_size = NULL){

  n = length(Z)
  H = length(methods.list.all)

  if(is.null(block.sum)){
    block.sum = summary_block(Z, block)
  }
  block = block.sum$block
  B = block.sum$B
  nb = block.sum$nb
  mb = block.sum$mb
  mb_ctrl = block.sum$mb_ctrl
  units.block = block.sum$units.block
  block.levels = block.sum$block.levels

  if(is.null(scores.list.all)){
    scores.list.all = list()
    for(h in 1:H){
      method.list.all = methods.list.all[[h]]
      scores.list.all[[h]] = score_all_blocks(nb, method.list.all)
    }
  }

  ## calculating mean and sigma of each stratified rank sum statistic
  mu_vec = rep(NA, H)
  sig_vec = rep(NA, H)
  for(h in 1:H){
    method.list.all = methods.list.all[[h]]
    score.list.all = scores.list.all[[h]]
    tmp = mu_sigma_single_weight_strat_rank_sum_stat(Z, block, method.list.all, score.list.all, weight, block.sum)

    mu_vec[h] = tmp$mu
    sig_vec[h] = tmp$sigma
  }

  ## generating random assignment vector corresponding to given block structure
  # Vectorized computation of test statistics
  # Key insight: for the null distribution we can take Y = 1:n, so within-block
  # ranks are fixed. Pre-compute weighted rank scores per unit, then use
  # matrix multiplication against permuted assignment matrices.

  # Pre-compute reordered scores for each method (n x H matrix)
  reordered_scores = matrix(0, nrow = n, ncol = H)
  for(h in 1:H){
    score.list.all = scores.list.all[[h]]
    for(i in 1:B){
      if(mb[i] == 0){
        next
      }
      units = units.block[[i]]
      reordered_scores[units, h] = weight[i] / mb[i] * score.list.all[[i]]
    }
  }

  if(!is.null(Z.perm)){
    test.stat.raw = crossprod(reordered_scores, Z.perm)
    test.stat.all = (test.stat.raw - mu_vec) / sig_vec
    return(apply(test.stat.all, 2, max))
  }

  if(is.null(chunk_size)){
    # Target ~5e6 entries per chunk for the n x chunk assignment matrix
    chunk_size = max(1L, min(null.max, as.integer(floor(5e6 / n))))
  }
  chunk_size = max(1L, min(null.max, as.integer(chunk_size)))

  combo_cache = new.env(parent = emptyenv())
  test.stat.max = numeric(null.max)
  for(offset in seq.int(1L, null.max, by = chunk_size)){
    size = min(chunk_size, null.max - offset + 1L)
    Z.chunk = .generate_Z_chunk_block(units.block = units.block, mb = mb, n = n, size = size, combo_cache = combo_cache)
    test.stat.raw = crossprod(reordered_scores, Z.chunk)
    test.stat.all = (test.stat.raw - mu_vec) / sig_vec
    test.stat.max[offset:(offset + size - 1L)] = apply(test.stat.all, 2, max)
  }

  return(test.stat.max)
}


#' Compute means and standard deviations for multiple statistics
#'
#' Calculates the mean and standard deviation for multiple weighted
#' stratified rank sum statistics.
#'
#' @param Z An n-dimensional binary treatment assignment vector.
#' @param block An n-dimensional vector specifying block membership.
#' @param weight A B-dimensional vector of block weights.
#' @param methods.list.all A list of lists of method specifications.
#' @param scores.list.all Optional pre-computed scores.
#' @param block.sum Optional pre-computed block summary.
#'
#' @return A list with components \code{mu} and \code{sigma}, each
#'   an H-dimensional vector for H statistics.
#'
#' @keywords internal
mu_sigma_list = function(Z, block, weight, methods.list.all, scores.list.all = NULL, block.sum = NULL){
  if(is.null(block.sum)){
    block.sum = summary_block(Z, block)
  }
  block = block.sum$block
  B = block.sum$B
  nb = block.sum$nb
  mb = block.sum$mb
  mb_ctrl = block.sum$mb_ctrl
  units.block = block.sum$units.block
  block.levels = block.sum$block.levels


  H = length(methods.list.all)
  mu = rep(0, H)
  sigma = rep(0, H)


  if(is.null(scores.list.all)){
    scores.list.all = list()
    for(h in 1:H){
      method.list.all = methods.list.all[[h]]
      scores.list.all[[h]] = score_all_blocks(nb, method.list.all)
    }
  }

  for(h in 1:H){
    method.list.all = methods.list.all[[h]]
    score.list.all = scores.list.all[[h]]
    temp = mu_sigma_single_weight_strat_rank_sum_stat(Z, block, method.list.all, score.list.all, weight, block.sum)
    mu[h] = temp$mu
    sigma[h] = temp$sigma
  }
  return(list(mu = mu, sigma = sigma))
}



#' Compute coefficient matrices for optimization
#'
#' Calculates the coefficient matrices needed for the optimization problem
#' in the combined stratified rank sum test.
#'
#' @param Z An n-dimensional binary treatment assignment vector.
#' @param Y An n-dimensional observed outcome vector.
#' @param block An n-dimensional vector specifying block membership.
#' @param c Threshold for the null hypothesis.
#' @param methods.list.all A list of lists of method specifications.
#' @param scores.list.all Optional pre-computed scores.
#' @param block.sum Optional pre-computed block summary.
#'
#' @return A list with H elements, each containing a list of B matrices.
#'   Each matrix is 2 x (mb+1) containing the number of units and minimum
#'   statistic values for different coverage scenarios.
#'
#' @keywords internal
comb_matrix_block <- function(Z, Y, block, c, methods.list.all, scores.list.all = NULL, block.sum = NULL){
  if(is.null(block.sum)){
    block.sum = summary_block(Z, block)
  }
  block = block.sum$block
  B = block.sum$B
  nb = block.sum$nb
  mb = block.sum$mb
  mb_ctrl = block.sum$mb_ctrl
  units.block = block.sum$units.block
  block.levels = block.sum$block.levels


  H = length(methods.list.all)

  N = length(Y)

  if(is.null(scores.list.all)){
    scores.list.all = list()
    for(h in 1:H){
      method.list.all = methods.list.all[[h]]
      scores.list.all[[h]] = score_all_blocks(nb, method.list.all)
    }
  }

  total.list = list()
  for (l in 1 : H){
    method.list.all = methods.list.all[[l]]
    score.list.all = scores.list.all[[l]]
    Tlist = list()
    for(i in 1:B){
      Zb = Z[ units.block[[i]] ]
      Yb = Y[ units.block[[i]] ]
      Ti = matrix(nrow = 2, ncol = mb[i] + 1)
      Ti[1,] = 0:mb[i]

      # The whole ii sequence comes from one ranking rather than from
      # mb[i] + 1 separate calls, each of which re-ran sort_treat and rank
      # over the same block.  See min_stat_path() in R/breakpoints.R.
      Ti[2, ] = min_stat_path(Zb, Yb, c, score.list.all[[i]])
      Tlist[[i]] = Ti
    }
    total.list[[l]] <- Tlist
  }
  return(total.list)
}


##############################################################################
## Stratum-level combination functions (comb.method = 2)
##############################################################################

#' Standard deviation for polynomial rank transformation
#'
#' Calculates the standard deviation of a polynomial rank score with
#' parameter r using the analytical formula.
#'
#' @param r Power parameter for polynomial rank transformation.
#'
#' @return Standard deviation value.
#'
#' @keywords internal
poly_sd <- function(r){
  sd = sqrt( (2*r-1)^(-1) - r^(-2) )
  return(sd)
}

#' Standard deviation of stratum-specific polynomial rank statistics
#'
#' Calculates the standard deviation for polynomial rank statistics
#' within each stratum.
#'
#' @param method.list.all A list of method specifications, one per block.
#'
#' @return A vector of standard deviations, one per stratum.
#'
#' @keywords internal
poly_strata_sd <- function(method.list.all){
  s = length(method.list.all)
  result = rep(NA, s)
  for(i in 1 : s){
    r = method.list.all[[i]]$r
    result[i] = poly_sd(r)
  }
  return(result)
}

#' Block-wise means and standard deviations for multiple methods
#'
#' Computes the mean and standard deviation for each block under
#' multiple rank score methods.
#'
#' @param Z An n-dimensional binary treatment assignment vector.
#' @param block An n-dimensional vector specifying block membership.
#' @param methods.list.all A list of H method specifications, each containing
#'   B block-specific methods.
#' @param scores.list.all Optional pre-computed scores.
#' @param weight A B-dimensional vector of block weights.
#' @param block.sum Optional pre-computed block summary.
#'
#' @return A list with:
#'   \item{mu_list}{List of B vectors, each of length H with block-wise means}
#'   \item{sig_list}{List of B vectors, each of length H with block-wise SDs}
#'
#' @keywords internal
mu_sd_block <- function(Z, block, methods.list.all, scores.list.all = NULL, weight,
                        block.sum = NULL){

  H = length(methods.list.all)

  if(is.null(block.sum)){
    block.sum = summary_block(Z, block)
  }
  block = block.sum$block
  B = block.sum$B
  nb = block.sum$nb
  mb = block.sum$mb
  mb_ctrl = block.sum$mb_ctrl
  units.block = block.sum$units.block
  block.levels = block.sum$block.levels


  if(is.null(scores.list.all)){
    scores.list.all = list()
    for(h in 1 : H){
      method.list.all = methods.list.all[[h]]
      scores.list.all[[h]] = score_all_blocks(nb, method.list.all)
    }
  }

  mu_list = list()
  sig_list = list()

  mu_tmp_mat = matrix(NA, nrow = H, ncol = B)
  sig_tmp_mat = matrix(NA, nrow = H, ncol = B)
  for(h in 1 : H){
    method.list.all = methods.list.all[[h]]
    score.list.all = scores.list.all[[h]]
    mu_tmp_mat[h , ] = mu_sigma_single_weight_strat_rank_sum_stat(Z, block, method.list.all, score.list.all, weight, block.sum)$mu_block
    sig_tmp_mat[h, ] = poly_strata_sd(method.list.all)
  }

  # arrange mean and standard deviations from multiple methods by blocks
  mu_list = list()
  sig_list = list()
  for(b in 1 : B){
    mu_list[[b]] = mu_tmp_mat[, b]
    sig_list[[b]] = sig_tmp_mat[, b]
  }
  return(list(mu_list = mu_list,
              sig_list = sig_list))
}

#' Compute stratum-wise standardized coefficient matrices
#'
#' Calculates the stratum-wise standardized test statistic values for
#' the stratum-level combination method (comb.method = 2).
#'
#' @param Z An n-dimensional binary treatment assignment vector.
#' @param Y An n-dimensional observed outcome vector.
#' @param block An n-dimensional vector specifying block membership.
#' @param c Threshold for the null hypothesis.
#' @param methods.list.all A list of H method specifications.
#' @param weight A B-dimensional vector of block weights.
#' @param scores.list.all Optional pre-computed scores.
#' @param block.sum Optional pre-computed block summary.
#'
#' @return A list with H elements, each containing B matrices with
#'   stratum-wise standardized test statistics.
#'
#' @keywords internal
comb_matrix_block_stratum <- function(Z, Y, block, c,
                                      methods.list.all,
                                      weight,
                                      scores.list.all = NULL, block.sum = NULL){

  H = length(methods.list.all)

  if(is.null(block.sum)){
    block.sum = summary_block(Z, block)
  }

  block = block.sum$block
  B = block.sum$B
  nb = block.sum$nb
  mb = block.sum$mb
  mb_ctrl = block.sum$mb_ctrl
  units.block = block.sum$units.block
  block.levels = block.sum$block.levels

  H = length(methods.list.all)
  N = length(Y)

  if(is.null(scores.list.all)){
    scores.list.all = list()
    for(h in 1:H){
      method.list.all = methods.list.all[[h]]
      scores.list.all[[h]] = score_all_blocks(nb, method.list.all)
    }
  }

  total.list = list()
  for (l in 1 : H){
    method.list.all = methods.list.all[[l]]
    score.list.all = scores.list.all[[l]]

    mu.vec = mu_sigma_single_weight_strat_rank_sum_stat(Z, block, method.list.all, weight = weight)$mu_block
    sd.vec = poly_strata_sd(method.list.all)

    Tlist = list()
    for(i in 1:B){
      Zb = Z[ units.block[[i]] ]
      Yb = Y[ units.block[[i]] ]
      # Columns are exempt counts, and the exempt count is the number of
      # treated units in the block whose effect exceeds c, so it stops at
      # mb[i].  min_stat exempts min(m, n - k) units (R/CMRSS_CRE.R:275),
      # so columns beyond mb[i] repeated the mb[i] column's value while
      # carrying a larger index into the solver's budget constraint.  They
      # were dominated rather than wrong, and dropping them changes no
      # answer; tests/testthat/test-comb-matrix-block-stratum-range.R
      # checks that across every budget.
      Ti = matrix(nrow = 2, ncol = mb[i] + 1)
      Ti[1,] = 0 : mb[i]

      score = score.list.all[[i]]

      # One ranking gives the whole column sequence; see min_stat_path().
      Ti[2,] = weight[i] *
        ( 1 / mb[i] * min_stat_path(Zb, Yb, c, score) - mu.vec[i] ) / sd.vec[i]
      Tlist[[i]] = Ti
    }
    total.list[[l]] <- Tlist
  }
  return(total.list)
}

#' Maximum of stratum-wise standardized statistics
#'
#' Computes the maximum of standardized statistics among each stratum
#' across multiple methods.
#'
#' @inheritParams comb_matrix_block_stratum
#'
#' @return A list of B matrices, one per block, containing the maximum
#'   standardized test statistic values.
#'
#' @keywords internal
max_comb_matrix_block_stratum <- function(Z, Y, block, c, methods.list.all, scores.list.all = NULL, block.sum = NULL,
                                          weight){

  # Forward the caller's precomputed scores and block summary rather than
  # passing NULL, which made comb_matrix_block_stratum rebuild both at every
  # threshold the confidence-interval search visits.
  tmp.list = comb_matrix_block_stratum(Z, Y, block, c, methods.list.all,
                                       scores.list.all = scores.list.all,
                                       block.sum = block.sum, weight = weight)
  s = length(tmp.list[[1]])  # number of blocks
  result.list = lapply(seq_len(s), function(j){
    mats = lapply(tmp.list, '[[', j)
    do.call(pmax, mats)
  })
  return(result.list)
}

#' Compute stratum-specific scores for multiple methods
#'
#' Generates score vectors for a single stratum under multiple methods.
#'
#' @param n Number of units in the stratum.
#' @param methods.list A list of H method specifications for this stratum.
#'
#' @return A list of H score vectors.
#'
#' @keywords internal
score_stratum <- function(n, methods.list){
  H = length(methods.list)
  score.list = list()
  for(h in 1 : H){
    method.list = methods.list[[h]]
    score.list[[h]] = rank_score(n, method.list)
  }
  return(score.list)
}

#' Stratum-specific weighted maximum rank sum statistic
#'
#' Calculates the weighted maximum rank sum statistic for a single stratum
#' across multiple methods.
#'
#' @param Z Treatment assignment vector for the stratum.
#' @param Y Outcome vector for the stratum.
#' @param methods.list A list of H method specifications.
#' @param score.list Optional pre-computed scores.
#' @param mu_vec Vector of H means.
#' @param sig_vec Vector of H standard deviations.
#' @param weight.b Weight for this stratum.
#'
#' @return Weighted maximum test statistic value.
#'
#' @keywords internal
single_weight_max_rank_sum_stat_stratum <- function(Z, Y,
                                                    methods.list,
                                                    score.list = NULL,
                                                    mu_vec,
                                                    sig_vec,
                                                    weight.b){
  H = length(methods.list)
  n = length(Z)
  m = sum(Z)

  if(is.null(score.list)){
    score.list = score_stratum(n, methods.list)
  }

  result_vec = rep(NA, H)
  for(h in 1 : H){
    score = score.list[[h]]
    mu = mu_vec[h]
    sd = sig_vec[h]
    rank_stat = ( 1 / m * sum( score[rank(Y, ties.method = "first")[Z == 1] ])  - mu ) / sd
    result_vec[h] = rank_stat
  }

  # taking maximum
  max_stat = max(result_vec)

  # adding weight
  result = weight.b * max_stat

  return(result)
}

#' Generate null distribution for stratum-level combination
#'
#' Generates the randomization null distribution for the stratum-level
#' combination method (comb.method = 2).
#'
#' @param Z An n-dimensional binary treatment assignment vector.
#' @param block An n-dimensional vector specifying block membership.
#' @param methods.list.all A list of H method specifications.
#' @param scores.list.all Optional pre-computed scores.
#' @param mu_sd_block_list Optional pre-computed means and SDs from mu_sd_block().
#' @param null.max Number of permutations.
#' @param weight A B-dimensional vector of block weights.
#' @param block.sum Optional pre-computed block summary.
#' @param Z.perm Optional pre-computed permutation matrix.
#' @param chunk_size Optional number of permutations per chunk when
#'   \code{Z.perm} is \code{NULL}. If \code{NULL}, chosen adaptively based on \code{n}.
#'
#' @return A numeric vector of the null distribution (length null.max).
#'
#' @keywords internal
com_null_dist_block_stratum <- function(Z, block, methods.list.all,
                                        scores.list.all = NULL,
                                        mu_sd_block_list = NULL,
                                        null.max = 10^4, weight,
                                        block.sum = NULL, Z.perm = NULL,
                                        chunk_size = NULL){

  H = length(methods.list.all)
  n = length(Z)

  if(is.null(block.sum)){
    block.sum = summary_block(Z, block)
  }

  block = block.sum$block
  B = block.sum$B
  nb = block.sum$nb
  mb = block.sum$mb
  mb_ctrl = block.sum$mb_ctrl
  units.block = block.sum$units.block
  block.levels = block.sum$block.levels

  if(is.null(scores.list.all)){
    scores.list.all = list()
    for(h in 1:H){
      method.list.all = methods.list.all[[h]]
      scores.list.all[[h]] = score_all_blocks(nb, method.list.all)
    }
  }

  if(is.null(mu_sd_block_list)){
    mu_sd_block_list = mu_sd_block(Z, block, methods.list.all, scores.list.all, weight, block.sum)
  }

  if(!is.null(Z.perm)){
    null.max = ncol(Z.perm)
  }

  stat.null = rep(NA, null.max)

  methods.list.by.block = list()
  for(b in 1:B){
    methods.list.b = list()
    for(h in 1:H){
      methods.list.b[[h]] = methods.list.all[[h]][[b]]
    }
    methods.list.by.block[[b]] = methods.list.b
  }

  score_mat.by.block = list()
  for(b in 1:B){
    score_mat.by.block[[b]] = do.call(cbind, lapply(1:H, function(h) scores.list.all[[h]][[b]]))
  }

  if(is.null(chunk_size)){
    chunk_size = max(1L, min(null.max, as.integer(floor(5e6 / n))))
  }
  chunk_size = max(1L, min(null.max, as.integer(chunk_size)))

  combo_cache = new.env(parent = emptyenv())
  for(offset in seq.int(1L, null.max, by = chunk_size)){
    size = min(chunk_size, null.max - offset + 1L)
    if(!is.null(Z.perm)){
      Z.chunk = Z.perm[, offset:(offset + size - 1L), drop = FALSE]
    } else {
      Z.chunk = .generate_Z_chunk_block(units.block = units.block, mb = mb, n = n, size = size, combo_cache = combo_cache)
    }

    result = numeric(size)
    for(b in 1:B){
      if(mb[b] == 0){
        next
      }
      units = units.block[[b]]
      Z.block = Z.chunk[units, , drop = FALSE]
      raw = crossprod(score_mat.by.block[[b]], Z.block) # H x size
      stat = raw / mb[b]
      stat = (stat - mu_sd_block_list$mu_list[[b]]) / mu_sd_block_list$sig_list[[b]]
      result = result + weight[b] * apply(stat, 2, max)
    }

    stat.null[offset:(offset + size - 1L)] = result
  }

  return(stat.null)
}


# Note: Solver functions (Gurobi_sol_com, HiGHS_sol_com, solve_optimization)
# are now defined in R/solvers.R to support multiple optimization backends.



#' Combined test statistic and combined p-value for randomization test
#' for quantiles of individual treatment effects (treated-only)
#'
#' Obtain the combined stratified rank sum statistic and combined p-value for
#' testing the treated-only null hypothesis
#' \eqn{H_{k,c}^{\mathrm{treat}}}: at least \eqn{k} of the \eqn{m = \sum_i Z_i}
#' treated units have individual treatment effect at least \eqn{c} (greater
#' alternative), or the analogous less / two-sided variants. The quantile
#' index \eqn{k} ranges over \eqn{1, \ldots, m}; values outside that range
#' make the underlying LP coverage constraint \eqn{p = m - k} infeasible and
#' raise an error.
#'
#' @param Z An \eqn{n} dimensional treatment assignment vector.
#' @param Y An \eqn{n} dimensional observed outcome vector.
#' @param k An integer between 1 and \eqn{\sum_i Z_i} (the number of treated
#'   units) specifying which quantile of the individual treatment effects
#'   among the treated is of interest.
#' @param c A numerical object specifying the threshold for the null
#'   hypothesis.
#' @param block An \eqn{n} dimensional vector specifying block of each
#'   units.
#' @param methods.list.all A list of lists of lists. Corresponds to the
#'   method for each stratum for each different stratified rank sum
#'   statistic.
#' @param weight.name Weighting method to be implemented. If
#'   "asymp.opt", asymptotically optimal scheme under a class of local
#'   alternatives is adjusted, if "dist.free", design-free scheme is
#'   adjusted.
#' @param stat.null An vector whose empirical distribution
#'   approximates the randomization distribution of the combined
#'   stratified rank sum statistic.
#' @param null.max A positive integer representing the number of
#'   permutations for approximating the randomization distribution of
#'   the rank sum statistic.
#' @param statistic logical; if TRUE (default), also prints the
#'   combined stratified rank sum statistic.
#' @param opt.method Optimization method. Options include:
#'   \itemize{
#'     \item "ILP_gurobi": Integer linear programming with Gurobi solver
#'     \item "LP_gurobi": Linear programming relaxation with Gurobi solver
#'     \item "ILP_highs": Integer linear programming with HiGHS solver (open-source)
#'     \item "LP_highs": Linear programming relaxation with HiGHS solver
#'     \item "ILP" or "ILP_auto": Integer LP with auto-selected solver
#'     \item "LP" or "LP_auto": Linear programming with auto-selected solver
#'   }
#'   HiGHS is recommended for users without a Gurobi license. Both solvers
#'   produce equivalent results.
#' @param comb.method Integer specifying the combination method:
#'   \itemize{
#'     \item 1: Combine statistics across strata first, then take maximum across methods (default)
#'     \item 2: Take maximum across methods within each stratum first, then combine across strata
#'   }
#'
#' @return The p-value (and test statistic) for testing the specified
#'   null hypothesis of interest.
#'
#' @examples
#' \dontrun{
#' # Load the electric teachers dataset
#' data(electric_teachers)
#'
#' # Set up treatment, outcome, and blocking variable
#' Z <- electric_teachers$TxAny
#' Y <- electric_teachers$gain
#' block <- factor(electric_teachers$Site)
#' n <- length(Y)
#' s <- length(levels(block))  # number of strata
#'
#' # Define polynomial rank statistics with different r values
#' # Each method must be specified for each stratum
#' r.vec <- c(2, 6, 10)
#' methods.list.all <- list()
#' for (j in seq_along(r.vec)) {
#'   methods.list.all[[j]] <- lapply(1:s, function(i) {
#'     list(name = "Polynomial", r = r.vec[j], std = TRUE, scale = FALSE)
#'   })
#' }
#'
#' # Test if the 90th percentile treated-only effect is <= 0
#' k <- floor(0.9 * sum(Z))
#' result <- pval_comb_block(Z, Y, k, c = 0, block, methods.list.all,
#'                           opt.method = "ILP", null.max = 10000)
#' print(result)
#'
#' # Using stratum-level combination method (comb.method = 2)
#' result2 <- pval_comb_block(Z, Y, k, c = 0, block, methods.list.all,
#'                            opt.method = "ILP", comb.method = 2)
#'
#' # Using LP relaxation (faster, conservative)
#' result3 <- pval_comb_block(Z, Y, k, c = 0, block, methods.list.all,
#'                            opt.method = "LP")
#' }
#'
#' @export
pval_comb_block <- function(Z, Y, k, c,
                            block,
                            methods.list.all,
                            weight.name = "asymp.opt",
                            stat.null = NULL,
                            null.max = 10^5,
                            statistic = TRUE,
                            opt.method = "ILP_auto",
                            comb.method = 1) {

  check_outcome_ties(Y)

  # Parse optimization method to get solver and exactness
  opt_params <- parse_opt_method(opt.method)
  exact <- opt_params$exact
  solver <- opt_params$solver

  H <- length(methods.list.all)
  N <- length(Y)
  m <- sum(Z)
  p <- m - k

  # Validate k against the current LP coverage constraint at line above.
  # With p = m - k, the LP is feasible only when k is in 1..m; outside this
  # range the LP is infeasible, the test statistic returns Inf, and the
  # function previously reported p.value = 0 silently (false rejection).
  # The message speaks to the user: it says what k counts and gives the
  # range, rather than citing a source line they cannot see.
  if (!is.numeric(k) || length(k) != 1L || k < 1 || k > m) {
    stop(sprintf(paste0(
      "k = %s, but pval_comb_block() tests the k-th smallest effect among ",
      "the treated units, so k must be between 1 and sum(Z) = %d. ",
      "For the 90th percentile of treated effects use k = floor(0.9 * sum(Z)). ",
      "For quantiles of all units' effects use ",
      "com_block_conf_quant_larger(set = \"all\")."),
      format(k), m
    ))
  }
  # tau_(2.5) is not a quantile, so a fractional k would answer a question
  # nobody can state; the LP would accept it without complaint.
  if (k != round(k)) {
    stop(sprintf("k = %s must be a whole number between 1 and sum(Z) = %d.",
                 format(k), m))
  }

  block.sum <- summary_block(Z, block)
  block <- block.sum$block

  weight <- weight_scheme(block.sum, weight.name)
  scores.list.all <- list()
  for (h in 1:H) {
    method.list.all <- methods.list.all[[h]]
    scores.list.all[[h]] <- score_all_blocks(block.sum$nb, method.list.all)
  }

  if (comb.method == 1) {
    # Method 1: Combine across strata first, then max across methods
    ms_list <- mu_sigma_list(Z, block, weight, methods.list.all, scores.list.all, block.sum)

    if (is.null(stat.null)) {
      stat.null <- com_null_dist_block(Z, block, methods.list.all, scores.list.all,
                                       null.max, weight, block.sum, Z.perm = NULL)
    }

    coeflists <- comb_matrix_block(Z, Y, block, c, methods.list.all, scores.list.all, block.sum)
    stat.min <- solve_optimization(Z, block, weight, coeflists, p, ms_list,
                                    exact, block.sum, solver)$obj
    pval <- mean(stat.null >= stat.min)

  } else {
    # Method 2: Max across methods within each stratum first, then combine
    mu_sd_block_list <- mu_sd_block(Z, block, methods.list.all,
                                    scores.list.all, weight, block.sum = NULL)

    if (is.null(stat.null)) {
      stat.null <- com_null_dist_block_stratum(Z, block, methods.list.all,
                                               scores.list.all, mu_sd_block_list,
                                               null.max, weight,
                                               block.sum, Z.perm = NULL)
    }

    coeflist <- max_comb_matrix_block_stratum(Z, Y, block, c, methods.list.all,
                                              scores.list.all, block.sum = NULL, weight = weight)
    stat.min <- solve_stratum_optimization(coeflist, p, exact = exact, solver = solver)$obj
    pval <- mean(stat.null >= stat.min)
  }

  if (statistic == TRUE) {
    result <- c(pval, stat.min)
    names(result) <- c("p.value", "test.stat")
    return(result)
  } else {
    return(pval)
  }
}



#' Helper function for Simultaneous Inference for multiple quantiles on SRE
#'
#' Output is a lower limit of prediction intervals for prespecified quantiles
#' under the all-units convention.  Used internally by
#' \code{\link{com_block_conf_quant_larger}}.
#'
#' @details
#' The hypothesis tested at each \code{k} is the all-units null
#' \eqn{H_{k, c}: \tau_{(k)} \le c}, where
#' \eqn{\tau_{(1)} \le \tau_{(2)} \le \ldots \le \tau_{(n)}} are the sorted
#' individual treatment effects across all \eqn{n} units.  The LP coverage
#' constraint sets \code{p <- n - k} (lines 1185 and 1248 of
#' \code{R/CMRSS_SRE.R}); \code{k} ranges over \eqn{1, \ldots, n} and the
#' loop fills \code{quantiles[k]} for \eqn{k = n - m, \ldots, n}, where
#' \eqn{m = \mathrm{sum}(Z)}, padding the remaining lower entries with
#' \code{quantiles[n - m]}.
#'
#' By the validity result extending \cite{ZL24quantile} (cf. main.tex
#' eq:H_kc_t and the discussion at lines 479 and 519-520), the inverted
#' \eqn{1 - \alpha} lower bound for \eqn{\tau_{(k)}} is simultaneously a
#' valid \eqn{1 - \alpha} lower prediction bound for
#' \eqn{\tau_{\mathrm{treat}(k - n_{\mathrm{control}})}}, where
#' \eqn{n_{\mathrm{control}} = n - m}.  Equivalently, the LP minimum at
#' all-units index \eqn{k_R} matches the LP minimum that
#' \code{\link{pval_comb_block}} (treated-only convention,
#' \code{p <- m - k}) would compute at \eqn{k_C = k_R - n_{\mathrm{control}}}.
#' The two functions therefore produce the same numerical confidence
#' bounds despite indexing their \code{k} differently.
#'
#' \code{\link{com_block_conf_quant_larger}} calls this function in three
#' modes: directly with \code{(Z, Y)} for \code{set = "treat"}, and with
#' \code{(1 - Z, -Y)} for \code{set = "control"} and the control half of
#' \code{set = "all"}.  Under the swap, individual effects are preserved via
#' the relabeling
#' \eqn{\tilde{\tau}_{\mathrm{treat}(k)} = \tau_{\mathrm{control}(k)}} (cf.
#' main.tex line 479), so the swap-space lower CIs returned by this function
#' are directly lower CIs on \eqn{\tau_{\mathrm{control}(k)}} in the original
#' problem -- no sign flip or index reversal is needed.
#'
#' @inheritParams pval_comb_block
#' @param k.vec Optional vector of specific quantile indices to compute.
#'   If NULL, computes for all quantiles.
#' @param tol Tolerance for root-finding algorithm.
#' @param alpha Significance level for confidence intervals.
#'
#' @return Vector of length \code{n} giving lower confidence limits for the
#'   all-units quantiles \eqn{\tau_{(1)}, \ldots, \tau_{(n)}}.  Only entries
#'   indexed by \eqn{k = n - m, \ldots, n} (where \eqn{m = \mathrm{sum}(Z)})
#'   are computed by the LP; lower entries are padded.
#'
#' @keywords internal
com_block_conf_quant_larger_trt <- function(Z, Y,
                                            block,
                                            k.vec = NULL,
                                            methods.list.all,
                                            weight.name = "asymp.opt",
                                            opt.method = "ILP_auto",
                                            comb.method = 1,
                                            stat.null = NULL, null.max = 10^4,
                                            Z.perm = NULL,
                                            tol = 0.01,
                                            alpha = 0.1) {

  # Parse optimization method to get solver and exactness
  opt_params <- parse_opt_method(opt.method)
  exact <- opt_params$exact
  solver <- opt_params$solver

  n <- length(Z)
  m <- sum(Z)

  block.sum <- summary_block(Z, block)
  block <- block.sum$block
  B <- block.sum$B

  H <- length(methods.list.all)

  weight <- weight_scheme(block.sum, weight.name)
  scores.list.all <- list()
  for (h in 1:H) {
    method.list.all <- methods.list.all[[h]]
    scores.list.all[[h]] <- score_all_blocks(block.sum$nb, method.list.all)
  }

  if (comb.method == 1) {
    ms_list <- mu_sigma_list(Z, block, weight, methods.list.all, scores.list.all, block.sum)

    if (is.null(stat.null)) {
      stat.null <- com_null_dist_block(Z, block, methods.list.all, scores.list.all,
                                       null.max, weight, block.sum, Z.perm = Z.perm)
    }
  } else {
    if (is.null(stat.null)) {
      stat.null <- com_null_dist_block_stratum(Z, block, methods.list.all,
                                               scores.list.all,
                                               mu_sd_block_list = NULL,
                                               null.max = null.max, weight = weight,
                                               block.sum = block.sum, Z.perm = Z.perm)
    }
  }

  # The critical value's position must count the simulated assignments
  # actually used. A caller-supplied Z.perm or stat.null can differ in length
  # from null.max; using null.max then picks the wrong order statistic (with
  # 2000 draws and null.max = 10^4, alpha = 0.10 became alpha = 0.50).
  null.max <- length(stat.null)
  thres <- sort(stat.null, decreasing = TRUE)[floor(null.max * alpha) + 1]

  Y1.max <- max(Y[Z == 1])
  Y1.min <- min(Y[Z == 1])
  Y0.max <- max(Y[Z == 0])
  Y0.min <- min(Y[Z == 0])

  c.max <- Y1.max - Y0.min + tol
  c.min <- Y1.min - Y0.max - tol

  if (is.null(k.vec)) {

    quantiles <- rep(NA, n)

    for (k in n:(n - m)) {

      if (k < n) {
        c.max <- quantiles[k + 1]
      }
      p <- n - k

      f <- function(c) {
        if (comb.method == 1) {
          coeflists <- comb_matrix_block(Z, Y, block, c, methods.list.all, scores.list.all, block.sum)
          stat.min <- solve_optimization(Z, block, weight, coeflists, p, ms_list,
                                          exact, block.sum, solver)$obj
        } else {
          coeflist <- max_comb_matrix_block_stratum(Z, Y, block, c, methods.list.all,
                                                    scores.list.all, block.sum, weight = weight)
          stat.min <- solve_stratum_optimization(coeflist, p, exact = exact, solver = solver)$obj
        }
        return(stat.min - thres)
      }

      tmp1 <- f(c.min)
      tmp2 <- f(c.max)
      if (tmp1 <= 0) {
        c.sol <- -Inf
      }
      if (tmp2 > 0) {
        c.sol <- c.max
      }
      if (tmp1 > 0 & tmp2 <= 0) {
        c.sol <- uniroot(f, interval = c(c.min, c.max), extendInt = "downX", tol = tol)$root
        c.sol <- round(c.sol, digits = -log10(tol))

        if (f(c.sol) <= 0) {
          while (f(c.sol) <= 0) {
            c.sol <- c.sol - tol
          }
          c.sol <- c.sol + tol
        } else {
          while (f(c.sol) > 0) {
            c.sol <- c.sol + tol
          }
        }
      }
      quantiles[k] <- c.sol

      if (c.sol == -Inf & k > (n - m)) {
        quantiles[(n - m):(k - 1)] <- -Inf
        break
      }
    }

    if (n - m > 1) {
      quantiles[1:(n - m - 1)] <- quantiles[n - m]
    }

    quantiles[quantiles > (Y1.max - Y0.min) + tol / 2] <- Inf

  } else {

    k.vec.sort <- sort(k.vec, decreasing = FALSE)
    j.max <- length(k.vec.sort)
    j.min <- max(sum(k.vec <= (n - m)), 1)

    quantiles <- rep(NA, j.max)

    for (j in j.max:j.min) {

      k <- k.vec.sort[j]
      p <- n - k
      if (j < j.max) {
        c.max <- quantiles[j + 1]
      }

      f <- function(c) {
        if (comb.method == 1) {
          coeflists <- comb_matrix_block(Z, Y, block, c, methods.list.all, scores.list.all, block.sum)
          stat.min <- solve_optimization(Z, block, weight, coeflists, p, ms_list,
                                          exact, block.sum, solver)$obj
        } else {
          coeflist <- max_comb_matrix_block_stratum(Z, Y, block, c, methods.list.all,
                                                    scores.list.all, block.sum = NULL, weight = weight)
          stat.min <- solve_stratum_optimization(coeflist, p, exact = exact, solver = solver)$obj
        }
        return(stat.min - thres)
      }

      tmp1 <- f(c.min)
      tmp2 <- f(c.max)
      if (tmp1 <= 0) {
        c.sol <- -Inf
      }
      if (tmp2 > 0) {
        c.sol <- c.max
      }
      if (f(c.min) > 0 & f(c.max) <= 0) {
        c.sol <- uniroot(f, interval = c(c.min, c.max), extendInt = "downX",
                         tol = tol)$root
        c.sol <- round(c.sol, digits = -log10(tol))
        if (f(c.sol) <= 0) {
          while (f(c.sol) <= 0) {
            c.sol <- c.sol - tol
          }
          c.sol <- c.sol + tol
        } else {
          while (f(c.sol) > 0) {
            c.sol <- c.sol + tol
          }
        }
      }
      quantiles[j] <- c.sol

      if (c.sol == -Inf & j > j.min) {
        quantiles[j.min:(j - 1)] <- -Inf
        break
      }
    }

    if (j.min > 1) {
      quantiles[1:(j.min - 1)] <- quantiles[j.min]
    }

  }
  return(quantiles)
}


#' Simultaneous bound for confidence interval using combined rank sum statistic
#' on stratified randomized experiment
#'
#' Computes simultaneous confidence intervals for quantiles of individual
#' treatment effects in stratified randomized experiments.
#'
#' @param Z An \eqn{n} dimensional treatment assignment vector.
#' @param Y An \eqn{n} dimensional observed outcome vector.
#' @param block An \eqn{n} dimensional vector specifying block of each
#'   unit.
#' @param set Set of quantiles of interest. Options:
#'   \itemize{
#'     \item "treat": Prediction intervals for effect quantiles among treated units
#'     \item "control": Prediction intervals for effect quantiles among control units
#'     \item "all": Confidence intervals for all effect quantiles
#'   }
#' @param methods.list.all A list of lists of lists. Corresponds to the
#'   method for each stratum for each different stratified rank sum
#'   statistic.
#' @param weight.name Weighting method to be implemented. If
#'   "asymp.opt", asymptotically optimal scheme under a class of local
#'   alternatives is adjusted, if "dist.free", design-free scheme is
#'   adjusted.
#' @param opt.method Optimization method. Options include:
#'   \itemize{
#'     \item "ILP_gurobi": Integer linear programming with Gurobi solver
#'     \item "LP_gurobi": Linear programming relaxation with Gurobi solver
#'     \item "ILP_highs": Integer linear programming with HiGHS solver (open-source)
#'     \item "LP_highs": Linear programming relaxation with HiGHS solver
#'     \item "ILP" or "ILP_auto": Integer LP with auto-selected solver
#'     \item "LP" or "LP_auto": Linear programming with auto-selected solver
#'   }
#'   HiGHS is recommended for users without a Gurobi license. Both solvers
#'   produce equivalent results.
#' @param comb.method Integer specifying the combination method:
#'   \itemize{
#'     \item 1: Combine statistics across strata first, then take maximum across methods (default)
#'     \item 2: Take maximum across methods within each stratum first, then combine across strata
#'   }
#' @param stat.null An vector whose empirical distribution
#'   approximates the randomization distribution of the combined
#'   stratified rank sum statistic. With \code{set = "control"} it must be
#'   the null distribution for the relabeled experiment, in which the
#'   control units are treated. It cannot be used with \code{set = "all"},
#'   because the two halves need different null distributions; supply
#'   \code{Z.perm} instead.
#' @param null.max A positive integer representing the number of
#'   permutations for approximating the randomization distribution of
#'   the rank sum statistic. Ignored when \code{Z.perm} or
#'   \code{stat.null} is supplied: the number of columns of \code{Z.perm},
#'   or the length of \code{stat.null}, is used instead.
#' @param Z.perm Optional pre-computed n x null.max matrix of permuted
#'   treatment assignments, in the original labeling (1 = treated). If
#'   provided, this matrix will be used instead of generating new
#'   permutations; for the control bounds the function uses
#'   \code{1 - Z.perm}. Can be generated using
#'   \code{assign_block(summary_block(Z, block), null.max)}.
#' @param tol A numerical object specifying the precision of the
#'   obtained confidence intervals. For example, if tol = 10^(-3),
#'   then the confidence limits are precise up to 3 digits.
#' @param alpha A numerical object, where 1-alpha indicates the
#'   confidence level.
#'
#' @return A vector specifying lower limits of prediction (confidence) intervals for
#'   quantiles k = 1 ~ m (or n - m, or n).
#'
#' @examples
#' \dontrun{
#' # Load the electric teachers dataset
#' data(electric_teachers)
#'
#' # Set up treatment, outcome, and blocking variable
#' Z <- electric_teachers$TxAny
#' Y <- electric_teachers$gain
#' block <- factor(electric_teachers$Site)
#' s <- length(levels(block))  # number of strata (7 sites)
#'
#' # Define polynomial rank statistics with different r values
#' # For SRE, methods must be specified for EACH stratum within EACH statistic
#' r.vec <- c(2, 6, 10)
#' methods.list.all <- list()
#' for (j in seq_along(r.vec)) {
#'   methods.list.all[[j]] <- lapply(1:s, function(i) {
#'     list(name = "Polynomial", r = r.vec[j], std = TRUE, scale = FALSE)
#'   })
#' }
#'
#' # Prediction intervals for effect quantiles among treated units
#' # Uses asymptotically optimal weighting across strata
#' ci.treat <- com_block_conf_quant_larger(Z, Y, block,
#'                                         set = "treat",
#'                                         methods.list.all = methods.list.all,
#'                                         weight.name = "asymp.opt",
#'                                         opt.method = "ILP",
#'                                         null.max = 10000,
#'                                         tol = 0.01,
#'                                         alpha = 0.05)
#'
#' # Prediction intervals for control units
#' ci.control <- com_block_conf_quant_larger(Z, Y, block,
#'                                           set = "control",
#'                                           methods.list.all = methods.list.all,
#'                                           opt.method = "ILP",
#'                                           alpha = 0.05)
#'
#' # Confidence intervals for all effect quantiles
#' ci.all <- com_block_conf_quant_larger(Z, Y, block,
#'                                       set = "all",
#'                                       methods.list.all = methods.list.all,
#'                                       opt.method = "ILP",
#'                                       alpha = 0.10)
#'
#' # Using stratum-level combination (comb.method = 2)
#' # Takes maximum across methods within each stratum first
#' ci.treat.m2 <- com_block_conf_quant_larger(Z, Y, block,
#'                                            set = "treat",
#'                                            methods.list.all = methods.list.all,
#'                                            opt.method = "ILP",
#'                                            comb.method = 2,
#'                                            alpha = 0.05)
#'
#' # Using LP relaxation for faster computation
#' ci.treat.lp <- com_block_conf_quant_larger(Z, Y, block,
#'                                            set = "treat",
#'                                            methods.list.all = methods.list.all,
#'                                            opt.method = "LP",
#'                                            alpha = 0.05)
#' }
#'
#' @seealso \code{\link{com_conf_quant_larger_cre}} for completely randomized experiments,
#'   \code{\link{pval_comb_block}} for hypothesis testing
#' @export
com_block_conf_quant_larger <- function(Z, Y,
                                        block,
                                        set = "all",
                                        methods.list.all = NULL,
                                        weight.name = "asymp.opt",
                                        opt.method = "ILP_auto",
                                        comb.method = 1,
                                        stat.null = NULL,
                                        null.max = 10^4,
                                        Z.perm = NULL,
                                        tol = 0.01,
                                        alpha = 0.1) {

  check_outcome_ties(Y)
  n <- length(Z)

  # The control bounds come from the relabeled experiment (Z <- 1 - Z,
  # Y <- -Y). Simulated assignments must be relabeled the same way, as
  # com_conf_quant_larger_cre does; otherwise every column scores the wrong
  # units.
  Z.perm.control <- if (is.null(Z.perm)) NULL else 1 - Z.perm

  # One null distribution cannot serve both halves: the relabeled experiment
  # has nb - mb treated units per stratum rather than mb.
  if (set == "all" && !is.null(stat.null)) {
    stop("stat.null cannot be used with set = \"all\": the treated and ",
         "control halves need different null distributions. Supply Z.perm ",
         "instead, or call set = \"treat\" and set = \"control\" separately.")
  }

  if (set == "treat") {
    ci.treat <- com_block_conf_quant_larger_trt(Z, Y, block,
                                                k.vec = NULL,
                                                methods.list.all,
                                                weight.name,
                                                opt.method,
                                                comb.method,
                                                stat.null, null.max, Z.perm, tol, alpha)
    ci.treat <- ci.treat[(n - sum(Z) + 1):n]
    return(ci.treat)
  }

  if (set == "control") {
    Z <- 1 - Z
    Y <- -Y
    ci.control <- com_block_conf_quant_larger_trt(Z, Y, block,
                                                  k.vec = NULL,
                                                  methods.list.all,
                                                  weight.name,
                                                  opt.method,
                                                  comb.method,
                                                  stat.null, null.max, Z.perm.control, tol,
                                                  alpha)
    ci.control <- ci.control[(n - sum(Z) + 1):n]
    return(ci.control)
  }

  if (set == "all") {
    ci.treat <- com_block_conf_quant_larger_trt(Z, Y, block,
                                                k.vec = NULL,
                                                methods.list.all,
                                                weight.name,
                                                opt.method,
                                                comb.method,
                                                stat.null, null.max, Z.perm, tol,
                                                alpha = alpha / 2)
    ci.treat <- ci.treat[(n - sum(Z) + 1):n]

    Z <- 1 - Z
    Y <- -Y

    ci.control <- com_block_conf_quant_larger_trt(Z, Y, block,
                                                  k.vec = NULL,
                                                  methods.list.all,
                                                  weight.name,
                                                  opt.method,
                                                  comb.method,
                                                  stat.null, null.max, Z.perm.control, tol,
                                                  alpha = alpha / 2)
    ci.control <- ci.control[(n - sum(Z) + 1):n]

    ci.all <- sort(c(ci.treat, ci.control))
    return(ci.all)
  }
}
