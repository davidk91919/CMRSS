# Null distribution for minimum p-values

Computes the distribution of the minimum p-value, particularly Monte
Carlo samples from the combined null distribution F in Theorem 1 of the
paper.

## Usage

``` r
comb_null_dist_cre(
  n,
  m,
  methods.list,
  Z.perm = NULL,
  nperm = 10^4,
  stat.null.mult = NULL
)
```

## Arguments

- n:

  Total number of units.

- m:

  Number of treated units.

- methods.list:

  A list of method specifications.

- Z.perm:

  Optional permutation matrix.

- nperm:

  Number of permutations for null distribution.

- stat.null.mult:

  Optional pre-computed null distribution matrix.

## Value

A numeric vector of minimum p-values under the null. Its length is the
number of draws actually available: `ncol(Z.perm)` when a permutation
matrix is supplied, `ncol(stat.null.mult)` when a null distribution is,
and `nperm` otherwise.

## Examples

``` r
set.seed(1)
methods.list <- list(list(name = "Wilcoxon", scale = FALSE),
                     list(name = "Stephenson", s = 3, scale = FALSE))
cnd <- comb_null_dist_cre(n = 20, m = 12, methods.list = methods.list,
                          nperm = 500)
quantile(cnd, c(0.05, 0.10))
#>     5%    10% 
#> 0.0459 0.0874 
```
