# Compute minimum test statistic

Calculates the minimum value of the rank sum statistic under the null
hypothesis, assuming at most n-k treated units can have effects greater
than c.

## Usage

``` r
min_stat(Z, Y, k, c, method.list = NULL, score = NULL, ind.sort.treat = NULL)
```

## Arguments

- Z:

  An n-dimensional binary treatment assignment vector.

- Y:

  An n-dimensional observed outcome vector.

- k:

  Quantile index (between 1 and n).

- c:

  Threshold for the null hypothesis.

- method.list:

  A list specifying the scoring method.

- score:

  Optional pre-computed score vector.

- ind.sort.treat:

  Optional pre-computed sorted treatment indices.

## Value

The minimum test statistic value.

## References

Caughey, D., Dafoe, A., Li, X., & Miratrix, L. (2023).

## Examples

``` r
Z <- c(0, 0, 0, 1, 1, 1)
Y <- c(1, 2, 3, 10, 20, 30)
score <- rank_score(6, list(name = "Wilcoxon", scale = FALSE))
# k = 6 exempts no treated unit, so this is the treated rank sum: 4 + 5 + 6
min_stat(Z, Y, k = 6, c = 0, score = score)
#> [1] 15
```
