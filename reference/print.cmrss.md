# Print a cmrss result

Reports how many of the sorted effects have a lower confidence bound
above `c`: if the bound for the k-th smallest of N effects is above c,
so are the bounds for the k-th through N-th, and at least N - k + 1
units have effects above c.

## Usage

``` r
# S3 method for class 'cmrss'
print(x, ...)
```

## Arguments

- x:

  An object returned by
  [`cmrss()`](https://bowers-illinois-edu.github.io/CMRSS/reference/cmrss.md).

- ...:

  Ignored.

## Value

`x`, invisibly.
