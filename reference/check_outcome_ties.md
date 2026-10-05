# Warn when ties in the outcome make results depend on row order

Warn when ties in the outcome make results depend on row order

## Usage

``` r
check_outcome_ties(Y, threshold = 0.05)
```

## Arguments

- Y:

  Outcome vector.

- threshold:

  Largest share of units that may share one value before the warning is
  given.

## Value

`Y`, invisibly. Called for its warning, which has class
`"cmrss_ties_warning"`.
