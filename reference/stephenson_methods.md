# Stephenson score lists, one per statistic and block

An s larger than a block's size is lowered to that size. The score
choose(r - 1, n_b - 1) is then 1 for the block's top rank and 0 below
it, so the block still counts instead of scoring 0 everywhere.

## Usage

``` r
stephenson_methods(s, nb)
```

## Arguments

- s:

  Stephenson parameters.

- nb:

  Block sizes (one value for a completely randomized experiment).

## Value

A list with one element per value of `s`, each a list with one score
specification per block, the form `methods.list.all` takes.
