# Default score parameters for cmrss(), polynomial zeta or Stephenson s

Three values: 2, s_max = min(floor(4 m / q_min), floor(n / 2)), and
their geometric middle round(sqrt(2 s_max)).

## Usage

``` r
default_score_parameters(n, m, alpha = 0.05)
```

## Arguments

- n:

  Number of units.

- m:

  Number of treated units.

- alpha:

  Level of the test.

## Value

A sorted integer vector of distinct values.
