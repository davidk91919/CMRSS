# Combining Multiple Rank Sum Statistics

> R package for conducting randomization inference for quantiles of
> individual treatment effects, using combined rank sum statistics, both
> for completely randomized and stratified randomized experiments

## Installation

``` r

devtools::install_github("jwbowers/CMRSS")
```

### Solver Dependencies

For stratified randomized experiments (SRE), the package requires an
optimization solver. You have two options:

**Option 1: HiGHS (Recommended - Open Source)**

``` r

install.packages("highs")
```

HiGHS is a high-performance open-source solver that requires no license.

**Option 2: Gurobi (Commercial)**

``` r

# Requires a Gurobi license (free for academics)
# See: https://www.gurobi.com/documentation/current/quickstart_mac/r_ins_the_r_package.html
```

Both solvers produce equivalent results. HiGHS is recommended for users
without a Gurobi license.

## Load the Package

``` r

library(CMRSS)
```

## Quick Start: Stratified Randomized Experiments

``` r

# Example: Testing treatment effects in a stratified experiment
set.seed(123)

# Create data
s <- 5  # number of strata
n <- 10 # units per stratum
N <- s * n

block <- factor(rep(1:s, each = n))
Z <- rep(0, N)
for (i in 1:s) {
  block_idx <- which(block == i)
  Z[sample(block_idx, n/2)] <- 1
}

Y0 <- rnorm(N)
Y1 <- Y0 + 1  # constant treatment effect of 1
Y <- Z * Y1 + (1 - Z) * Y0

# Set up methods (Wilcoxon scores)
methods.list.all <- list()
methods.list.all[[1]] <- lapply(1:s, function(i) list(name = "Wilcoxon", scale = FALSE))

# Compute p-value using HiGHS solver
result <- pval_comb_block(Z, Y, k = floor(0.9 * N), c = 0,
                          block, methods.list.all,
                          opt.method = "ILP_highs")
print(result)
```

## Solver Options

The `opt.method` parameter controls which optimization solver to use:

| opt.method | Description |
|----|----|
| `"ILP_highs"` | Integer linear programming with HiGHS (open-source) |
| `"LP_highs"` | Linear programming relaxation with HiGHS |
| `"ILP_gurobi"` | Integer linear programming with Gurobi |
| `"LP_gurobi"` | Linear programming relaxation with Gurobi |
| `"ILP"` or `"ILP_auto"` | Auto-select available solver (prefers Gurobi) |
| `"LP"` or `"LP_auto"` | Auto-select for LP relaxation |

### Checking Solver Availability

``` r

# Check which solvers are installed
solver_available("highs")   # TRUE if highs package is installed
solver_available("gurobi")  # TRUE if gurobi package is installed

# Get the default solver
get_default_solver()  # Returns "highs" or "gurobi"
```

### Comparing Solvers

Both HiGHS and Gurobi produce equivalent results. You can verify this:

``` r

# Same results with different solvers
result_highs <- pval_comb_block(Z, Y, k, c, block, methods.list.all,
                                opt.method = "ILP_highs")

result_gurobi <- pval_comb_block(Z, Y, k, c, block, methods.list.all,
                                 opt.method = "ILP_gurobi")

# These should be equal (up to numerical precision)
all.equal(result_highs, result_gurobi)
```

## Main Functions

### Completely Randomized Experiments (CRE)

``` r

?comb_p_val_cre              # Get p-value from combined statistics in CRE
?com_conf_quant_larger_cre   # Get confidence intervals for quantiles in CRE
```

### Stratified Randomized Experiments (SRE)

``` r

?pval_comb_block             # Get p-value from combined statistics in SRE
?com_block_conf_quant_larger # Get confidence intervals for quantiles in SRE
```

### Solver Utilities

``` r

?solve_optimization   # Unified optimization interface
?solver_available     # Check if a solver is available
?get_default_solver   # Get the default solver
?parse_opt_method     # Parse optimization method string
```

## Examples

### Example 1: P-value Calculation

``` r

library(CMRSS)

# Setup
set.seed(42)
s <- 4; n <- 8; N <- s * n
block <- factor(rep(1:s, each = n))
Z <- unlist(lapply(1:s, function(i) sample(c(rep(1, n/2), rep(0, n/2)))))
Y <- rnorm(N) + Z * 0.8  # Treatment effect of 0.8

# Methods: combine Wilcoxon and Stephenson scores
methods.list.all <- list(
  lapply(1:s, function(i) list(name = "Wilcoxon", scale = FALSE)),
  lapply(1:s, function(i) list(name = "Stephenson", s = 3, scale = FALSE))
)

# Test H0: tau_(k) <= 0 at 90th percentile
result <- pval_comb_block(Z, Y, k = floor(0.9 * N), c = 0,
                          block, methods.list.all,
                          opt.method = "ILP_highs")
cat("P-value:", result["p.value"], "\n")
cat("Test statistic:", result["test.stat"], "\n")
```

### Example 2: Confidence Intervals

``` r

# Compute simultaneous confidence intervals for treatment effect quantiles
ci <- com_block_conf_quant_larger(Z, Y, block,
                                  set = "treat",  # for treated units
                                  methods.list.all = methods.list.all,
                                  opt.method = "ILP_highs",
                                  alpha = 0.1)
print(ci)
```

## References

- Kim, D. and Li, X. (2024). Combining Multiple Rank Sum Statistics for
  Inference on Treatment Effect Quantiles.

## License

MIT License
