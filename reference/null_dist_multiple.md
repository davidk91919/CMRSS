# Generate null distributions for multiple rank sum statistics

Computes the randomization null distribution for multiple rank sum
statistics simultaneously under complete randomization.

## Usage

``` r
null_dist_multiple(
  n,
  m,
  methods.list = NULL,
  nperm = 10^5,
  Z.perm = NULL,
  chunk_size = 1000
)
```

## Arguments

- n:

  Total number of units.

- m:

  Number of treated units.

- methods.list:

  A list of method specifications.

- nperm:

  Number of permutations.

- Z.perm:

  Optional pre-computed permutation matrix.

- chunk_size:

  Integer chunk size used to generate permutations when `Z.perm` is
  `NULL`. Smaller values use less memory.

## Value

An H x nperm matrix where H is the number of statistics.

## Examples

``` r
set.seed(1)
methods.list <- list(list(name = "Wilcoxon", scale = FALSE),
                     list(name = "Stephenson", s = 3, scale = FALSE))
dim(null_dist_multiple(n = 20, m = 12, methods.list = methods.list,
                       nperm = 500))
#> [1]   2 500
```
