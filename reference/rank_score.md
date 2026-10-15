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

  - `divide_by_top`: (for Stephenson) logical. TRUE, and also when left
    out, divides every score by the top unit's score choose(n - 1, s -
    1), so scores run from 0 to 1. FALSE gives the binomial coefficients
    choose(r - 1, s - 1) themselves.

  - `r`: (for Polynomial) the power parameter

  - `std`: (for Polynomial) logical, use Puri(1965) normalization

  - `scale`: logical, standardize scores to mean 0 and sd 1

## Value

A numeric vector of length n containing the rank scores.

## Details

Polynomial scores optionally use Puri(1965) normalization. Stephenson
scores use binomial coefficients, by default divided by the top unit's
score. The division is computed on the log scale, exp(lchoose(r - 1,
s - 1) - lchoose(n - 1, s - 1)), because the coefficients themselves
exceed the largest double once they pass about 1e308: with s = 500, from
n = 1,032 on, where the undivided scores are Inf and every p-value is
NA. Without strata the division changes no p-value, since the observed
statistic and every null draw are multiplied by the same constant. With
strata of unequal size it changes the test, because each stratum's
scores are divided by a different number: undivided, a stratum of 100
with s = 10 has top score choose(99, 9) = 1.7e12 and a stratum of 20 has
choose(19, 9) = 92,378, so the larger stratum carries nearly all the
weight. Wilcoxon scores are simply the ranks 1 to n.

## Examples

``` r
rank_score(5, list(name = "Wilcoxon", scale = FALSE))
#> [1] 1 2 3 4 5
rank_score(5, list(name = "Stephenson", s = 3, scale = FALSE))
#> [1] 0.0000000 0.0000000 0.1666667 0.5000000 1.0000000
rank_score(5, list(name = "Stephenson", s = 3, divide_by_top = FALSE, scale = FALSE))
#> [1] 0 0 1 3 6
rank_score(5, list(name = "Polynomial", r = 2, std = TRUE, scale = FALSE))
#> [1] 0.1666667 0.3333333 0.5000000 0.6666667 0.8333333
```
