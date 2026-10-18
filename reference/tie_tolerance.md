# Margin within which two numbers count as tied

Margin within which two numbers count as tied

## Usage

``` r
tie_tolerance(x)
```

## Arguments

- x:

  The observed number (or numbers) the comparison is made against.

## Value

`1e-9 * pmax(1, abs(x))`.
