# Calculate rank scores

Computes rank scores for rank-based test statistics. Supports Wilcoxon,
Stephenson, and Polynomial rank scores.

## Usage

``` r
rank_score(
  n,
  method.list = list(name = "Polynomial", r, std = TRUE, scale = FALSE)
)
```

## Arguments

- n:

  Number of units.

- method.list:

  A list specifying the scoring method:

  - `name`: "Wilcoxon", "Stephenson", or "Polynomial"

  - `s`: (for Stephenson) the parameter s

  - `r`: (for Polynomial) the power parameter

  - `std`: (for Polynomial) logical, use Puri(1965) normalization

  - `scale`: logical, standardize scores to mean 0 and sd 1

## Value

A numeric vector of length n containing the rank scores.

## Details

Polynomial scores optionally use Puri(1965) normalization. Stephenson
scores use binomial coefficients. Wilcoxon scores are simply the ranks 1
to n.

## Examples

``` r
rank_score(5, list(name = "Wilcoxon", scale = FALSE))
#> [1] 1 2 3 4 5
rank_score(5, list(name = "Stephenson", s = 3, scale = FALSE))
#> [1] 0 0 1 3 6
rank_score(5, list(name = "Polynomial", r = 2, std = TRUE, scale = FALSE))
#> [1] 0.1666667 0.3333333 0.5000000 0.6666667 0.8333333
```
