# Generate null distribution for rank sum statistic

Computes the randomization null distribution of a single rank sum
statistic under complete randomization.

## Usage

``` r
null_dist(
  n,
  m,
  method.list = NULL,
  score = NULL,
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

- method.list:

  A list specifying the scoring method.

- score:

  Optional pre-computed score vector.

- nperm:

  Number of permutations.

- Z.perm:

  Optional pre-computed permutation matrix.

- chunk_size:

  Integer chunk size used to generate permutations when `Z.perm` is
  `NULL`. Smaller values use less memory.

## Value

A numeric vector of length nperm containing the null distribution.

## Examples

``` r
set.seed(1)
score <- rank_score(20, list(name = "Wilcoxon", scale = FALSE))
nd <- null_dist(n = 20, m = 12, score = score, nperm = 500)
quantile(nd, c(0.5, 0.95))
#>    50%    95% 
#> 127.00 148.05 
```
