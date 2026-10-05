# Polynomial score lists, one per statistic and block

A polynomial score (r / (n_b + 1))^(zeta - 1) is never 0, so zeta is
used unchanged in every block.

## Usage

``` r
polynomial_methods(zeta, nb)
```

## Arguments

- zeta:

  Polynomial parameters.

- nb:

  Block sizes (one value for a completely randomized experiment).

## Value

A list in the form `methods.list.all` takes.
