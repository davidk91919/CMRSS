# Share of draws at or below a value, counting near-ties as ties

Share of draws at or below a value, counting near-ties as ties

## Usage

``` r
share_at_most(draws, value)
```

## Arguments

- draws:

  Numeric vector of null draws.

- value:

  A single observed value.

## Value

The share of `draws` that are at most `value + tie_tolerance(value)`.
