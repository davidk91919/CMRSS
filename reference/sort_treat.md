# Sort treated units by outcome rank

Returns the indices of treated units sorted by their outcome values in
increasing order.

## Usage

``` r
sort_treat(Y, Z)
```

## Arguments

- Y:

  An n-dimensional observed outcome vector.

- Z:

  An n-dimensional binary treatment assignment vector.

## Value

Integer vector of indices of treated units, sorted by increasing
outcome.

## Examples

``` r
Z <- c(0, 0, 0, 1, 1, 1)
Y <- c(1, 2, 3, 30, 10, 20)
sort_treat(Y, Z)   # 5, 6, 4: the treated units ordered by outcome
#> [1] 5 6 4
```
