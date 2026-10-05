# Simultaneous bound for confidence interval using combined rank sum statistic on stratified randomized experiment

Computes simultaneous confidence intervals for quantiles of individual
treatment effects in stratified randomized experiments.

## Usage

``` r
com_block_conf_quant_larger(
  Z,
  Y,
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
  alpha = 0.1
)
```

## Arguments

- Z:

  An \\n\\ dimensional treatment assignment vector.

- Y:

  An \\n\\ dimensional observed outcome vector.

- block:

  An \\n\\ dimensional vector specifying block of each unit.

- set:

  Set of quantiles of interest. Options:

  - "treat": Prediction intervals for effect quantiles among treated
    units

  - "control": Prediction intervals for effect quantiles among control
    units

  - "all": Confidence intervals for all effect quantiles

- methods.list.all:

  A list of lists of lists. Corresponds to the method for each stratum
  for each different stratified rank sum statistic.

- weight.name:

  Weighting method to be implemented. If "asymp.opt", asymptotically
  optimal scheme under a class of local alternatives is adjusted, if
  "dist.free", design-free scheme is adjusted.

- opt.method:

  Optimization method. Options include:

  - "ILP_gurobi": Integer linear programming with Gurobi solver

  - "LP_gurobi": Linear programming relaxation with Gurobi solver

  - "ILP_highs": Integer linear programming with HiGHS solver
    (open-source)

  - "LP_highs": Linear programming relaxation with HiGHS solver

  - "ILP" or "ILP_auto": Integer LP with auto-selected solver

  - "LP" or "LP_auto": Linear programming with auto-selected solver

  HiGHS is recommended for users without a Gurobi license. Both solvers
  produce equivalent results.

- comb.method:

  Integer specifying the combination method:

  - 1: Combine statistics across strata first, then take maximum across
    methods (default)

  - 2: Take maximum across methods within each stratum first, then
    combine across strata

- stat.null:

  An vector whose empirical distribution approximates the randomization
  distribution of the combined stratified rank sum statistic. With
  `set = "control"` it must be the null distribution for the relabeled
  experiment, in which the control units are treated. It cannot be used
  with `set = "all"`, because the two halves need different null
  distributions; supply `Z.perm` instead.

- null.max:

  A positive integer representing the number of permutations for
  approximating the randomization distribution of the rank sum
  statistic. Ignored when `Z.perm` or `stat.null` is supplied: the
  number of columns of `Z.perm`, or the length of `stat.null`, is used
  instead.

- Z.perm:

  Optional pre-computed n x null.max matrix of permuted treatment
  assignments, in the original labeling (1 = treated). If provided, this
  matrix will be used instead of generating new permutations; for the
  control bounds the function uses `1 - Z.perm`. Can be generated using
  `assign_block(summary_block(Z, block), null.max)`.

- tol:

  A numerical object specifying the precision of the obtained confidence
  intervals. For example, if tol = 10^(-3), then the confidence limits
  are precise up to 3 digits.

- alpha:

  A numerical object, where 1-alpha indicates the confidence level.

## Value

A vector specifying lower limits of prediction (confidence) intervals
for quantiles k = 1 ~ m (or n - m, or n).

## See also

[`com_conf_quant_larger_cre`](https://bowers-illinois-edu.github.io/CMRSS/reference/com_conf_quant_larger_cre.md)
for completely randomized experiments,
[`pval_comb_block`](https://bowers-illinois-edu.github.io/CMRSS/reference/pval_comb_block.md)
for hypothesis testing

## Examples

``` r
if (FALSE) { # \dontrun{
# Load the electric teachers dataset
data(electric_teachers)

# Set up treatment, outcome, and blocking variable
Z <- electric_teachers$TxAny
Y <- electric_teachers$gain
block <- factor(electric_teachers$Site)
s <- length(levels(block))  # number of strata (7 sites)

# Define polynomial rank statistics with different r values
# For SRE, methods must be specified for EACH stratum within EACH statistic
r.vec <- c(2, 6, 10)
methods.list.all <- list()
for (j in seq_along(r.vec)) {
  methods.list.all[[j]] <- lapply(1:s, function(i) {
    list(name = "Polynomial", r = r.vec[j], std = TRUE, scale = FALSE)
  })
}

# Prediction intervals for effect quantiles among treated units
# Uses asymptotically optimal weighting across strata
ci.treat <- com_block_conf_quant_larger(Z, Y, block,
                                        set = "treat",
                                        methods.list.all = methods.list.all,
                                        weight.name = "asymp.opt",
                                        opt.method = "ILP",
                                        null.max = 10000,
                                        tol = 0.01,
                                        alpha = 0.05)

# Prediction intervals for control units
ci.control <- com_block_conf_quant_larger(Z, Y, block,
                                          set = "control",
                                          methods.list.all = methods.list.all,
                                          opt.method = "ILP",
                                          alpha = 0.05)

# Confidence intervals for all effect quantiles
ci.all <- com_block_conf_quant_larger(Z, Y, block,
                                      set = "all",
                                      methods.list.all = methods.list.all,
                                      opt.method = "ILP",
                                      alpha = 0.10)

# Using stratum-level combination (comb.method = 2)
# Takes maximum across methods within each stratum first
ci.treat.m2 <- com_block_conf_quant_larger(Z, Y, block,
                                           set = "treat",
                                           methods.list.all = methods.list.all,
                                           opt.method = "ILP",
                                           comb.method = 2,
                                           alpha = 0.05)

# Using LP relaxation for faster computation
ci.treat.lp <- com_block_conf_quant_larger(Z, Y, block,
                                           set = "treat",
                                           methods.list.all = methods.list.all,
                                           opt.method = "LP",
                                           alpha = 0.05)
} # }
```
