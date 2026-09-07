# Combined test statistic and combined p-value for randomization test for quantiles of individual treatment effects (treated-only)

Obtain the combined stratified rank sum statistic and combined p-value
for testing the treated-only null hypothesis
\\H\_{k,c}^{\mathrm{treat}}\\: at least \\k\\ of the \\m = \sum_i Z_i\\
treated units have individual treatment effect at least \\c\\ (greater
alternative), or the analogous less / two-sided variants. The quantile
index \\k\\ ranges over \\1, \ldots, m\\; values outside that range make
the underlying LP coverage constraint \\p = m - k\\ infeasible and raise
an error.

## Usage

``` r
pval_comb_block(
  Z,
  Y,
  k,
  c,
  block,
  methods.list.all,
  weight.name = "asymp.opt",
  stat.null = NULL,
  null.max = 10^5,
  statistic = TRUE,
  opt.method = "ILP_auto",
  comb.method = 1
)
```

## Arguments

- Z:

  An \\n\\ dimensional treatment assignment vector.

- Y:

  An \\n\\ dimensional observed outcome vector.

- k:

  An integer between 1 and \\\sum_i Z_i\\ (the number of treated units)
  specifying which quantile of the individual treatment effects among
  the treated is of interest.

- c:

  A numerical object specifying the threshold for the null hypothesis.

- block:

  An \\n\\ dimensional vector specifying block of each units.

- methods.list.all:

  A list of lists of lists. Corresponds to the method for each stratum
  for each different stratified rank sum statistic.

- weight.name:

  Weighting method to be implemented. If "asymp.opt", asymptotically
  optimal scheme under a class of local alternatives is adjusted, if
  "dist.free", design-free scheme is adjusted.

- stat.null:

  An vector whose empirical distribution approximates the randomization
  distribution of the combined stratified rank sum statistic.

- null.max:

  A positive integer representing the number of permutations for
  approximating the randomization distribution of the rank sum
  statistic.

- statistic:

  logical; if TRUE (default), also prints the combined stratified rank
  sum statistic.

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

## Value

The p-value (and test statistic) for testing the specified null
hypothesis of interest.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load the electric teachers dataset
data(electric_teachers)

# Set up treatment, outcome, and blocking variable
Z <- electric_teachers$TxAny
Y <- electric_teachers$gain
block <- factor(electric_teachers$Site)
n <- length(Y)
s <- length(levels(block))  # number of strata

# Define polynomial rank statistics with different r values
# Each method must be specified for each stratum
r.vec <- c(2, 6, 10)
methods.list.all <- list()
for (j in seq_along(r.vec)) {
  methods.list.all[[j]] <- lapply(1:s, function(i) {
    list(name = "Polynomial", r = r.vec[j], std = TRUE, scale = FALSE)
  })
}

# Test if the 90th percentile treated-only effect is <= 0
k <- floor(0.9 * sum(Z))
result <- pval_comb_block(Z, Y, k, c = 0, block, methods.list.all,
                          opt.method = "ILP", null.max = 10000)
print(result)

# Using stratum-level combination method (comb.method = 2)
result2 <- pval_comb_block(Z, Y, k, c = 0, block, methods.list.all,
                           opt.method = "ILP", comb.method = 2)

# Using LP relaxation (faster, conservative)
result3 <- pval_comb_block(Z, Y, k, c = 0, block, methods.list.all,
                           opt.method = "LP")
} # }
```
