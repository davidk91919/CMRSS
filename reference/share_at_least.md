# Share of draws at or above a value, counting near-ties as ties

Share of draws at or above a value, counting near-ties as ties

## Usage

``` r
share_at_least(draws, value)
```

## Arguments

- draws:

  Numeric vector of null draws.

- value:

  A single observed value.

## Value

The share of `draws` that are at least `value - tie_tolerance(value)`.
