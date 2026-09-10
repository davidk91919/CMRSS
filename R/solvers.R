#' Solver Utilities for CMRSS Package
#'
#' This file contains optimization solver functions for solving the integer/linear
#' programming problems in stratified randomized experiments. It supports both
#' Gurobi and HiGHS solvers.
#'
#' @description
#' The package supports two optimization solvers:
#' \itemize{
#'   \item \strong{Gurobi}: Commercial solver (requires license). Fast and robust.
#'   \item \strong{HiGHS}: Open-source solver. No license required, good performance.
#' }
#'
#' Both solvers produce equivalent results for the same optimization problem.
#' HiGHS is recommended for users without a Gurobi license.

#' Check if a solver is available
#'
#' @param solver Character string: "gurobi" or "highs"
#' @return Logical indicating if the solver is available
#' @examples
#' solver_available("highs")
#' solver_available("gurobi")
#' @export
solver_available <- function(solver) {
  solver <- tolower(solver)
  if (solver == "gurobi") {
    return(requireNamespace("gurobi", quietly = TRUE))
  } else if (solver == "highs") {
    return(requireNamespace("highs", quietly = TRUE))
  } else {
    stop("Unknown solver: ", solver, ". Must be 'gurobi' or 'highs'.")
  }
}

#' Get default solver
#'
#' Returns the default solver based on availability. Prefers Gurobi if available
#' (to maintain backward compatibility with earlier versions/workflows), and
#' falls back to HiGHS.
#'
#' @return Character string indicating the available solver
#' @examples
#' get_default_solver()
#' @export
get_default_solver <- function() {
  if (solver_available("gurobi")) {
    return("gurobi")
  } else if (solver_available("highs")) {
    return("highs")
  } else {
    stop("No solver available. Please install either 'highs' (recommended, open-source) or 'gurobi'.\n",
         "To install highs: install.packages('highs')\n",
         "For gurobi, see: https://www.gurobi.com/documentation/current/quickstart_mac/r_ins_the_r_package.html")
  }
}


#' Optimization function using HiGHS solver
#'
#' Solves the integer/linear programming problem for computing the minimum test
#' statistic in stratified randomized experiments using the HiGHS solver.
#'
#' @param Z An n-dimensional treatment assignment vector.
#' @param block An n-dimensional vector specifying block of each unit.
#' @param weight A B-dimensional vector of block weights.
#' @param coeflists A list of H elements, each containing B matrices with
#'   block-specific test statistics for different values of k.
#' @param p Upper bound on the number of units with effects greater than c.
#' @param ms_list A list containing mu (means) and sigma (standard deviations)
#'   for each test statistic.
#' @param exact Logical; if TRUE, solve as integer linear program (ILP).
#'   If FALSE, solve as linear program (LP relaxation).
#' @param block.sum Optional pre-computed block summary from summary_block().
#'
#' @return A list with:
#'   \item{sol}{Solution vector}
#'   \item{obj}{Optimal objective value}
#'
#' @details
#' This function formulates and solves an optimization problem to find the
#' minimum value of the combined test statistic under the constraint that
#' at most p units have effects greater than a threshold c.
#'
#' The problem has the following structure:
#' \itemize{
#'   \item Variables: x_sj (binary/continuous), eta_h, theta
#'   \item Objective: minimize theta
#'   \item Constraints: probability constraints per stratum, coverage constraint,
#'     test statistic equality constraints, and normalization constraints
#' }
#'
#' @seealso \code{\link{Gurobi_sol_com}} for the Gurobi implementation,
#'   \code{\link{solve_optimization}} for the unified wrapper
#'
#' @examples
#' \dontrun{
#' # This is typically called internally by pval_comb_block()
#' # Direct usage requires pre-computed coefficient lists
#' result <- HiGHS_sol_com(Z, block, weight, coeflists, p, ms_list, exact = TRUE)
#' }
#' @keywords internal
#' @importFrom Matrix sparseMatrix
HiGHS_sol_com <- function(Z, block, weight, coeflists, p, ms_list, exact = TRUE, block.sum = NULL) {

  if (!requireNamespace("highs", quietly = TRUE)) {
    stop("Package 'highs' is required but not installed.\n",
         "Install it with: install.packages('highs')")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required but not installed.\n",
         "Install it with: install.packages('Matrix')")
  }

  ## nn: number of all units + B
  ## B: number of strata
  ## H: number of stat to combine
  ## Variables are organized by: x's (nn), eta's (H), theta (1)
  ## weight is a B-dimensional vector

  H <- length(coeflists)
  B <- length(coeflists[[1]])

  if (is.null(block.sum)) {
    block.sum <- summary_block(Z, block)
  }
  block <- block.sum$block
  mb <- block.sum$mb
  mb_ctrl <- block.sum$mb_ctrl
  units.block <- block.sum$units.block
  block.levels <- block.sum$block.levels

  mmb <- mb + 1  # 0 ~ m_s + 1; for the constraint of sum of the number of k's in each strata

  ## Total number of x variables
  nn <- sum(mmb)

  ## Total number of variables: x's (nn) + eta's (H) + theta (1)
  n_vars <- nn + H + 1

  ## Objective: minimize theta (last variable)
  obj <- rep(0, n_vars)
  obj[n_vars] <- 1

  ## Build constraint matrix using sparse representation
  ## Constraint types:
  ## 1. Sum of x in each strata = 1 (B constraints)
  ## 2. Sum of k*x <= p (1 constraint)
  ## 3. Test statistic constraints (H constraints)
  ## 4. Final normalization constraints (H constraints)

  Ai_parts <- list()
  Aj_parts <- list()
  Ax_parts <- list()

  ## 1. Sum of x in each strata = 1 (B equality constraints)
  Ai_parts[[length(Ai_parts) + 1L]] <- rep(1:B, mmb)
  Aj_parts[[length(Aj_parts) + 1L]] <- 1:nn
  Ax_parts[[length(Ax_parts) + 1L]] <- rep(1, nn)

  ## 2. Coverage constraint: sum of k*x <= p
  Ai_parts[[length(Ai_parts) + 1L]] <- rep(B + 1, nn)
  Aj_parts[[length(Aj_parts) + 1L]] <- 1:nn
  coeflist <- coeflists[[1]]
  k_coefs <- unlist(lapply(coeflist, function(mat) mat[1, ]), use.names = FALSE)
  Ax_parts[[length(Ax_parts) + 1L]] <- k_coefs

  ## 3. Test statistic constraints: T_h - eta_h = 0 (H constraints)
  for (h in 1:H) {
    coeflist <- coeflists[[h]]
    Ai_parts[[length(Ai_parts) + 1L]] <- rep(B + 1 + h, nn)
    Aj_parts[[length(Aj_parts) + 1L]] <- 1:nn
    t_coefs <- unlist(lapply(1:B, function(b) weight[b] * (1 / mb[b]) * coeflist[[b]][2, ]),
                      use.names = FALSE)
    Ax_parts[[length(Ax_parts) + 1L]] <- t_coefs
  }

  ## Add -1 coefficient for eta_h in each test statistic constraint
  Ai_parts[[length(Ai_parts) + 1L]] <- (B + 2):(B + 1 + H)
  Aj_parts[[length(Aj_parts) + 1L]] <- (nn + 1):(nn + H)
  Ax_parts[[length(Ax_parts) + 1L]] <- rep(-1, H)

  ## 4. Normalization constraints: theta - (1/sigma_h)*mu_h - eta_h <= 0 (H constraints)
  mu <- ms_list$mu
  sig_rev <- 1 / ms_list$sigma

  ## theta coefficient (-1) in each normalization constraint
  Ai_parts[[length(Ai_parts) + 1L]] <- (B + H + 2):(B + 2 * H + 1)
  Aj_parts[[length(Aj_parts) + 1L]] <- rep(n_vars, H)
  Ax_parts[[length(Ax_parts) + 1L]] <- rep(-1, H)

  ## eta_h coefficient (1/sigma_h) in each normalization constraint
  Ai_parts[[length(Ai_parts) + 1L]] <- (B + H + 2):(B + 2 * H + 1)
  Aj_parts[[length(Aj_parts) + 1L]] <- (nn + 1):(nn + H)
  Ax_parts[[length(Ax_parts) + 1L]] <- sig_rev

  Ai <- unlist(Ai_parts, use.names = FALSE)
  Aj <- unlist(Aj_parts, use.names = FALSE)
  Ax <- unlist(Ax_parts, use.names = FALSE)

  ## Create sparse constraint matrix
  n_constraints <- B + 1 + H + H
  A <- Matrix::sparseMatrix(i = Ai, j = Aj, x = Ax, dims = c(n_constraints, n_vars))

  ## Right-hand side values
  ## Constraints 1 (B equalities): rhs = 1
  ## Constraint 2 (coverage): rhs = p
  ## Constraints 3 (H equalities): rhs = 0
  ## Constraints 4 (H inequalities): rhs = mu * sig_rev (i.e., mu/sigma)
  rhs <- c(rep(1, B), p, rep(0, H), mu * sig_rev)

  ## Constraint directions for HiGHS: lhs <= Ax <= rhs
  ## For equality: lhs = rhs
  ## For <=: lhs = -Inf
  lhs <- c(rep(1, B), -Inf, rep(0, H), rep(-Inf, H))

  ## Variable bounds
  ## x's: [0, 1] for LP, binary for ILP
  ## eta's: [-Inf, Inf]
  ## theta: [-Inf, Inf]
  lower_bounds <- c(rep(0, nn), rep(-Inf, H + 1))
  upper_bounds <- c(rep(1, nn), rep(Inf, H + 1))

  ## Variable types
  ## HiGHS uses: "C" = continuous, "I" = integer
  ## Binary variables are integers with 0-1 bounds
  if (exact) {
    types <- c(rep("I", nn), rep("C", H + 1))
  } else {
    types <- rep("C", n_vars)
  }

  ## Solve with HiGHS
  result <- highs::highs_solve(
    L = obj,
    lower = lower_bounds,
    upper = upper_bounds,
    A = A,
    lhs = lhs,
    rhs = rhs,
    types = types,
    maximum = FALSE,  # minimize
    control = highs::highs_control(log_to_console = FALSE)
  )

  ## Check if solution was found
  if (result$status != 7) {  # 7 = kOptimal in HiGHS
    warning("HiGHS solver did not find optimal solution. Status: ", result$status_message)
  }

  return(list(sol = result$primal_solution, obj = result$objective_value))
}


#' Optimization function using Gurobi solver (original implementation)
#'
#' Solves the integer/linear programming problem for computing the minimum test
#' statistic in stratified randomized experiments using the Gurobi solver.
#'
#' @inheritParams HiGHS_sol_com
#'
#' @return A list with:
#'   \item{sol}{Solution vector}
#'   \item{obj}{Optimal objective value}
#'
#' @details
#' This is the original Gurobi-based implementation. For users without a Gurobi
#' license, \code{\link{HiGHS_sol_com}} provides equivalent functionality using
#' the open-source HiGHS solver.
#'
#' @seealso \code{\link{HiGHS_sol_com}} for the HiGHS implementation,
#'   \code{\link{solve_optimization}} for the unified wrapper
#'
#' @examples
#' \dontrun{
#' # This is typically called internally by pval_comb_block()
#' result <- Gurobi_sol_com(Z, block, weight, coeflists, p, ms_list, exact = TRUE)
#' }
#' @keywords internal
Gurobi_sol_com <- function(Z, block, weight, coeflists, p, ms_list, exact = TRUE, block.sum = NULL) {

  if (!requireNamespace("gurobi", quietly = TRUE)) {
    stop("Package 'gurobi' is required but not installed.\n",
         "See: https://www.gurobi.com/documentation/current/quickstart_mac/r_ins_the_r_package.html")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required but not installed.\n",
         "Install it with: install.packages('Matrix')")
  }

  ## nn: number of all units + B
  ## B: number of strata
  ## H: number of stat to combine
  ## Q is organized by x's(n*H), eta's(H), theta
  ## weight is a H by B matrix, where each hth row is a weight vector for hth statistic.

  model <- list()
  H <- length(coeflists)
  B <- length(coeflists[[1]])

  if (is.null(block.sum)) {
    block.sum <- summary_block(Z, block)
  }
  block <- block.sum$block
  mb <- block.sum$mb
  mb_ctrl <- block.sum$mb_ctrl
  units.block <- block.sum$units.block
  block.levels <- block.sum$block.levels

  mmb <- mb + 1  ## 0 ~ m_s + 1; for the constraint of sum of the number of k's in each strata


  ## parameter setting in the gurobi function

  nn <- sum(mmb)
  Q <- rep(0, nn + H + 1)
  Q[nn + H + 1] <- 1

  ## sum of x in each strata is 1
  Ai_parts <- list(rep(1:B, mmb))
  Aj_parts <- list(1:nn)
  x_parts <- list(rep(1, nn))

  # cost is p for all choices of stat and profit for h's stat is eta'h

  ## choosing x_{sj} for each stratum s, where x_{sj}=1 if and only if j units in this stratum have effects greater than c
  Ai_parts[[length(Ai_parts) + 1L]] <- rep((B + 1), nn)
  Aj_parts[[length(Aj_parts) + 1L]] <- 1:nn
  coeflist <- coeflists[[1]]
  x_parts[[length(x_parts) + 1L]] <- unlist(lapply(coeflist, function(mat) mat[1, ]), use.names = FALSE)

  ## listing t_s^h(j)s
  for (h in 1:H) {
    coeflist <- coeflists[[h]]
    Ai_parts[[length(Ai_parts) + 1L]] <- rep(B + 1 + h, nn)
    Aj_parts[[length(Aj_parts) + 1L]] <- 1:nn
    x_parts[[length(x_parts) + 1L]] <- unlist(lapply(1:B, function(b) weight[b] * (1 / mb[b]) * coeflist[[b]][2, ]),
                                       use.names = FALSE)
  }

  Ai_parts[[length(Ai_parts) + 1L]] <- (B + 1 + 1):(B + 1 + H)
  Aj_parts[[length(Aj_parts) + 1L]] <- (nn + 1):(nn + H)
  x_parts[[length(x_parts) + 1L]] <- rep(-1, H)

  mu <- ms_list$mu
  sig_rev <- 1 / ms_list$sigma

  ## for final constraints

  Ai_parts[[length(Ai_parts) + 1L]] <- rep((B + H + 1 + 1):(B + 2 * H + 1), 2)
  Aj_parts[[length(Aj_parts) + 1L]] <- c(rep(nn + H + 1, H), (nn + 1):(nn + H))
  x_parts[[length(x_parts) + 1L]] <- c(rep(-1, H), sig_rev)

  Ai <- unlist(Ai_parts, use.names = FALSE)
  Aj <- unlist(Aj_parts, use.names = FALSE)
  x <- unlist(x_parts, use.names = FALSE)

  A <- Matrix::sparseMatrix(Ai, Aj, x = x)


  model$A <- A
  model$obj <- Q
  model$modelsense <- "min"
  model$rhs <- c(rep(1, B), p, rep(0, H), mu * sig_rev)
  model$sense <- c(rep("=", B), "<=", rep("=", H), rep("<=", H))
  model$lb <- c(rep(0, nn), rep(-Inf, H + 1))
  if (exact) {
    model$vtype <- c(rep("B", nn), rep("C", H + 1))
  }


  params <- add_thread_param(list(OutputFlag = 0))
  result <- gurobi::gurobi(model, params)

  return(list(sol = result$x, obj = result$objval))
}


#' Add a thread cap to a Gurobi parameter list
#'
#' Gurobi's default is to use every core on the machine for one solve. That is
#' the right default for a caller solving one problem at a time, and the wrong
#' one for a caller that is already running in parallel: a `foreach` loop over
#' 14 workers, each solving with 14 threads, puts roughly 150 runnable threads
#' on 14 cores and the workers queue behind each other instead of solving.
#'
#' Only the caller knows whether it is already parallel, so the thread count is
#' the caller's to set, through `options(CMRSS.gurobi.threads = 1)`. Forked
#' workers inherit R options, so setting it once before a parallel loop reaches
#' every worker. Unset, which is the default, this returns `params` untouched
#' and Gurobi chooses as it always has.
#'
#' Threads is a resource setting rather than a modelling one, so capping it
#' does not change the solution; `tests/testthat/test-gurobi-threads.R` checks
#' that on a four-block problem.
#'
#' @param params A list of Gurobi parameters.
#'
#' @return `params`, with `Threads` added when the option is set.
#'
#' @keywords internal
add_thread_param <- function(params) {
  threads <- getOption("CMRSS.gurobi.threads", NULL)
  if (is.null(threads)) {
    return(params)
  }
  if (!is.numeric(threads) || length(threads) != 1L || !is.finite(threads) ||
      threads < 1) {
    stop("option CMRSS.gurobi.threads must be a single positive integer, not ",
         deparse(threads), call. = FALSE)
  }
  params$Threads <- as.integer(threads)
  params
}


#' Unified solver wrapper
#'
#' A unified interface for solving the optimization problem using either
#' Gurobi or HiGHS solver.
#'
#' @inheritParams HiGHS_sol_com
#' @param solver Character string specifying the solver to use. Options are:
#'   \itemize{
#'     \item "auto": Automatically select available solver (default)
#'     \item "highs": Use HiGHS solver (open-source)
#'     \item "gurobi": Use Gurobi solver (commercial, requires license)
#'   }
#'
#' @return A list with:
#'   \item{sol}{Solution vector}
#'   \item{obj}{Optimal objective value}
#'   \item{solver}{The solver that was used}
#'
#' @details
#' This function provides a unified interface for both solvers. When solver = "auto",
#' it will prefer Gurobi if available (for backward compatibility), otherwise fall back
#' to HiGHS.
#'
#' Both solvers produce equivalent results for the same optimization problem,
#' though there may be minor numerical differences (typically < 1e-6) due to

#' different internal algorithms.
#'
#' @seealso \code{\link{HiGHS_sol_com}}, \code{\link{Gurobi_sol_com}},
#'   \code{\link{solver_available}}, \code{\link{get_default_solver}}
#'
#' @examples
#' \dontrun{
#' # Automatically select solver
#' result <- solve_optimization(Z, block, weight, coeflists, p, ms_list)
#'
#' # Explicitly use HiGHS
#' result <- solve_optimization(Z, block, weight, coeflists, p, ms_list, solver = "highs")
#'
#' # Explicitly use Gurobi
#' result <- solve_optimization(Z, block, weight, coeflists, p, ms_list, solver = "gurobi")
#' }
#' @export
solve_optimization <- function(Z, block, weight, coeflists, p, ms_list,
                                exact = TRUE, block.sum = NULL, solver = "auto") {

  ## Determine which solver to use
  if (solver == "auto") {
    solver <- get_default_solver()
  } else {
    solver <- tolower(solver)
    if (!solver %in% c("gurobi", "highs")) {
      stop("Unknown solver: ", solver, ". Must be 'auto', 'gurobi', or 'highs'.")
    }
    if (!solver_available(solver)) {
      stop("Solver '", solver, "' is not available. Please install the corresponding package.")
    }
  }

  ## Call the appropriate solver
  if (solver == "highs") {
    result <- HiGHS_sol_com(Z, block, weight, coeflists, p, ms_list, exact, block.sum)
  } else {
    result <- Gurobi_sol_com(Z, block, weight, coeflists, p, ms_list, exact, block.sum)
  }

  result$solver <- solver
  return(result)
}


#' Parse optimization method string
#'
#' Parses the opt.method parameter to determine solver and exactness.
#'
#' @param opt.method Character string specifying optimization method:
#'   \itemize{
#'     \item "ILP_gurobi": Integer linear programming with Gurobi
#'     \item "LP_gurobi": Linear programming relaxation with Gurobi
#'     \item "ILP_highs": Integer linear programming with HiGHS
#'     \item "LP_highs": Linear programming relaxation with HiGHS
#'     \item "ILP_auto" or "ILP": Integer linear programming with auto-selected solver
#'     \item "LP_auto" or "LP": Linear programming with auto-selected solver
#'   }
#'
#' @return A list with:
#'   \item{exact}{Logical indicating if ILP (TRUE) or LP (FALSE)}
#'   \item{solver}{Character string: "gurobi", "highs", or "auto"}
#'
#' @examples
#' parse_opt_method("ILP_gurobi")  # list(exact = TRUE, solver = "gurobi")
#' parse_opt_method("LP_highs")    # list(exact = FALSE, solver = "highs")
#' parse_opt_method("ILP")         # list(exact = TRUE, solver = "auto")
#' @export
parse_opt_method <- function(opt.method) {
  opt.method <- toupper(opt.method)

  # Determine if ILP or LP
  if (grepl("^ILP", opt.method)) {
    exact <- TRUE
  } else if (grepl("^LP", opt.method)) {
    exact <- FALSE
  } else {
    stop("Invalid opt.method: ", opt.method,
         ". Must start with 'ILP' or 'LP'.")
  }

  # Determine solver
  if (grepl("GUROBI", opt.method)) {
    solver <- "gurobi"
  } else if (grepl("HIGHS", opt.method)) {
    solver <- "highs"
  } else if (grepl("AUTO", opt.method) || opt.method %in% c("ILP", "LP")) {
    solver <- "auto"
  } else {
    # Default: treat unknown suffix as auto
    solver <- "auto"
  }

  return(list(exact = exact, solver = solver))
}


##############################################################################
## Stratum-level combination solvers (for comb.method = 2)
##############################################################################

#' Stratum-level optimization using Gurobi solver
#'
#' Solves the optimization problem for the stratum-level combination method
#' (comb.method = 2) using Gurobi.
#'
#' @param coeflist A list of B matrices containing the stratum-wise
#'   standardized test statistics for different coverage scenarios.
#' @param p Upper bound on the number of units with effects greater than c.
#' @param exact Logical; if TRUE, solve as integer linear program (ILP).
#'   If FALSE, solve as linear program (LP relaxation).
#'
#' @return A list with:
#'   \item{sol}{Solution vector}
#'   \item{obj}{Optimal objective value}
#'
#' @seealso \code{\link{HiGHS_sol_stratum_com}} for the HiGHS implementation
#' @keywords internal
Gurobi_sol_stratum_com <- function(coeflist, p, exact = TRUE) {

  if (!requireNamespace("gurobi", quietly = TRUE)) {
    stop("Package 'gurobi' is required but not installed.\n",
         "See: https://www.gurobi.com/documentation/current/quickstart_mac/r_ins_the_r_package.html")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required but not installed.\n",
         "Install it with: install.packages('Matrix')")
  }

  model <- list()
  B <- length(coeflist)
  nb <- vapply(coeflist, ncol, integer(1))
  Q <- unlist(lapply(coeflist, function(mat) mat[2, ]), use.names = FALSE)

  n <- length(Q)
  indx <- c(0, cumsum(nb))
  Ai <- rep(B + 1, 2 * n)
  Aj <- rep(1:n, 2)
  x <- rep(1, 2 * n)
  for (i in 1:B) {
    Ai[(indx[i] + 1):indx[i + 1]] <- i
    x[(n + indx[i] + 1):(n + indx[i + 1])] <- coeflist[[i]][1, ]
  }
  A <- Matrix::sparseMatrix(Ai, Aj, x = x)
  model$A <- A
  model$obj <- Q
  model$modelsense <- "min"
  model$rhs <- c(rep(1, B), p)
  model$sense <- c(rep("=", B), "<=")
  if (exact) {
    model$vtype <- "B"
  }
  params <- add_thread_param(list(OutputFlag = 0))
  result <- gurobi::gurobi(model, params)
  return(list(sol = result$x, obj = result$objval))
}


#' Stratum-level optimization using HiGHS solver
#'
#' Solves the optimization problem for the stratum-level combination method
#' (comb.method = 2) using HiGHS.
#'
#' @inheritParams Gurobi_sol_stratum_com
#'
#' @return A list with:
#'   \item{sol}{Solution vector}
#'   \item{obj}{Optimal objective value}
#'
#' @seealso \code{\link{Gurobi_sol_stratum_com}} for the Gurobi implementation
#'
#' @keywords internal
HiGHS_sol_stratum_com <- function(coeflist, p, exact = TRUE) {

  if (!requireNamespace("highs", quietly = TRUE)) {
    stop("Package 'highs' is required but not installed.\n",
         "Install it with: install.packages('highs')")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required but not installed.\n",
         "Install it with: install.packages('Matrix')")
  }

  B <- length(coeflist)
  nb <- vector(length = B)

  for (i in 1:B) {
    nb[i] <- ncol(coeflist[[i]])
  }
  # Build objective vector from test statistics
  Q <- unlist(lapply(coeflist, function(mat) mat[2, ]), use.names = FALSE)

  n <- length(Q)

  # Build constraint matrix
  indx <- c(0, cumsum(nb))
  Ai <- rep(B + 1, 2 * n)
  Aj <- rep(1:n, 2)
  Ax <- rep(1, 2 * n)
  for (i in 1:B) {
    Ai[(indx[i] + 1):indx[i + 1]] <- i
    Ax[(n + indx[i] + 1):(n + indx[i + 1])] <- coeflist[[i]][1, ]
  }
  A <- Matrix::sparseMatrix(i = Ai, j = Aj, x = Ax, dims = c(B + 1, n))

  # Right-hand side and bounds
  rhs <- c(rep(1, B), p)
  lhs <- c(rep(1, B), -Inf)

  # Variable bounds and types
  lower_bounds <- rep(0, n)
  upper_bounds <- rep(1, n)

  if (exact) {
    types <- rep("I", n)
  } else {
    types <- rep("C", n)
  }

  # Solve with HiGHS
  result <- highs::highs_solve(
    L = Q,
    lower = lower_bounds,
    upper = upper_bounds,
    A = A,
    lhs = lhs,
    rhs = rhs,
    types = types,
    maximum = FALSE,
    control = highs::highs_control(log_to_console = FALSE)
  )

  if (result$status != 7) {
    warning("HiGHS solver did not find optimal solution. Status: ", result$status_message)
  }

  return(list(sol = result$primal_solution, obj = result$objective_value))
}


#' Unified solver wrapper for stratum-level optimization
#'
#' A unified interface for solving the stratum-level optimization problem
#' using either Gurobi or HiGHS solver.
#'
#' @inheritParams Gurobi_sol_stratum_com
#' @param solver Character string specifying the solver to use: "auto",
#'   "highs", or "gurobi".
#'
#' @return A list with:
#'   \item{sol}{Solution vector}
#'   \item{obj}{Optimal objective value}
#'   \item{solver}{The solver that was used}
#'
#' @seealso \code{\link{solve_optimization}} for the standard optimization
#'
#' @keywords internal
solve_stratum_optimization <- function(coeflist, p, exact = TRUE, solver = "auto") {

  # Determine which solver to use
  if (solver == "auto") {
    solver <- get_default_solver()
  } else {
    solver <- tolower(solver)
    if (!solver %in% c("gurobi", "highs")) {
      stop("Unknown solver: ", solver, ". Must be 'auto', 'gurobi', or 'highs'.")
    }
    if (!solver_available(solver)) {
      stop("Solver '", solver, "' is not available. Please install the corresponding package.")
    }
  }

  # Call the appropriate solver
  if (solver == "highs") {
    result <- HiGHS_sol_stratum_com(coeflist, p, exact)
  } else {
    result <- Gurobi_sol_stratum_com(coeflist, p, exact)
  }

  result$solver <- solver
  return(result)
}
