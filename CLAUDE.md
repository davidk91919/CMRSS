# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## Package Overview

CMRSS (Combining Multiple Rank Sum Statistics) is an R package for
randomization inference about quantiles of individual treatment effects.
It implements methods from Kim and Li (2024) for combining multiple rank
sum statistics to improve power when testing hypotheses about treatment
effect quantiles.

The package supports: - **Completely Randomized Experiments (CRE)**:
Functions prefixed with `_cre` - **Stratified Randomized Experiments
(SRE)**: Functions prefixed with `_block`

## Build and Test Commands

``` bash
# Run all tests
R CMD check .

# Run tests with testthat (faster iteration)
Rscript -e "devtools::test()"

# Run a specific test file
Rscript -e "testthat::test_file('tests/testthat/test-solvers.R')"

# Load package for development (uses devtools)
Rscript load.R

# Generate documentation
Rscript -e "roxygen2::roxygenize()"

# Install from local source
R CMD INSTALL .
```

## Architecture

### Core Module Structure

    R/
    ├── CMRSS_CRE.R      # Completely randomized experiments
    ├── CMRSS_SRE.R      # Stratified randomized experiments (main file)
    ├── solvers.R        # Optimization solver abstraction (HiGHS/Gurobi)
    ├── comparison_methods.R  # Benchmark methods from literature (RIQITE wrappers)
    ├── data.R           # Dataset documentation (electric_teachers)
    └── CMRSS-package.R  # Package-level documentation

### Key Functions

**CRE (Completely Randomized Experiments)**: -
[`comb_p_val_cre()`](https://bowers-illinois-edu.github.io/CMRSS/reference/comb_p_val_cre.md) -
Combined p-value from multiple rank statistics -
[`com_conf_quant_larger_cre()`](https://bowers-illinois-edu.github.io/CMRSS/reference/com_conf_quant_larger_cre.md) -
Confidence intervals for effect quantiles

**SRE (Stratified Randomized Experiments)**: -
[`pval_comb_block()`](https://bowers-illinois-edu.github.io/CMRSS/reference/pval_comb_block.md) -
Combined p-value with stratification -
[`com_block_conf_quant_larger()`](https://bowers-illinois-edu.github.io/CMRSS/reference/com_block_conf_quant_larger.md) -
Confidence intervals with stratification

**Solver Functions (R/solvers.R)**: -
[`solve_optimization()`](https://bowers-illinois-edu.github.io/CMRSS/reference/solve_optimization.md) -
Unified solver interface (auto-selects HiGHS or Gurobi) -
[`solver_available()`](https://bowers-illinois-edu.github.io/CMRSS/reference/solver_available.md) -
Check if a solver is installed -
[`parse_opt_method()`](https://bowers-illinois-edu.github.io/CMRSS/reference/parse_opt_method.md) -
Parse optimization method strings like “ILP_highs”

### Optimization Solvers

The package requires an optimization solver for SRE methods: - **HiGHS**
(recommended, open-source): `install.packages("highs")` - **Gurobi**
(commercial, requires license)

The `opt.method` parameter controls solver selection: - `"ILP_highs"` /
`"LP_highs"` - Use HiGHS - `"ILP_gurobi"` / `"LP_gurobi"` - Use Gurobi -
`"ILP"` / `"LP"` or `"ILP_auto"` / `"LP_auto"` - Auto-select

### Rank Score Methods

Methods are specified as lists with a `name` field: -
`list(name = "Wilcoxon", scale = FALSE)` - Standard Wilcoxon scores -
`list(name = "Stephenson", s = 3, scale = FALSE)` - Stephenson scores
with parameter s -
`list(name = "Polynomial", r = 2, std = TRUE, scale = FALSE)` -
Polynomial scores

For SRE, methods must be specified per stratum in a nested list
structure:

``` r

# H statistics, each with B block-specific methods
methods.list.all <- list(
  lapply(1:B, function(i) list(name = "Polynomial", r = 2, std = TRUE, scale = FALSE)),
  lapply(1:B, function(i) list(name = "Polynomial", r = 6, std = TRUE, scale = FALSE))
)
```

### Dependencies

- Required: stats, utils, Matrix, RIQITE (from GitHub:
  li-xinran/RIQITE), extraDistr
- Suggested: testthat, highs, gurobi, randomizr, RItools

## Testing Patterns

Tests use testthat edition 3. Test files cover: - `test-solvers.R` -
Solver availability and equivalence - `test-CMRSS_CRE.R` /
`test-CMRSS_SRE.R` - Main function tests - `test-rank_score.R` - Rank
score computation - `test-comparison_methods.R` - Benchmark method
wrappers

Example test pattern:

``` r

test_that("description", {
  skip_if_not(solver_available("highs"), "HiGHS not available")
  # test code
})
```

## Data

The package includes `electric_teachers` dataset - a stratified
randomized experiment from Heller et al. (2010) with 233 teachers across
7 sites.
