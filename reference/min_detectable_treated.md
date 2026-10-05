# Fewest treated units whose effects could be detected at level alpha

The smallest q for which the probability that the q largest outcomes all
belong to treated units, under no effect, is at most `alpha`: (m/n)((m -
1)/(n - 1)) ... ((m - q + 1)/(n - q + 1)) \<= alpha.

## Usage

``` r
min_detectable_treated(n, m, alpha = 0.05)
```

## Arguments

- n:

  Number of units.

- m:

  Number of treated units.

- alpha:

  Level of the test.

## Value

An integer, or `NA` when no q up to m reaches `alpha`.
