# Every exempt count's statistic from a single ranking

Returns `min_stat` for a block at exempt counts \\0, 1, \ldots, m_b\\,
the whole sequence `comb_matrix_block` needs, computed from one call to
`rank` instead of \\m_b + 1\\ of them.

## Usage

``` r
min_stat_path(Zb, Yb, c, score)
```

## Arguments

- Zb:

  Treatment assignment within one block.

- Yb:

  Observed outcomes within the same block.

- c:

  The threshold.

- score:

  Rank score vector of length `length(Yb)`.

## Value

A numeric vector of length `sum(Zb) + 1` giving the statistic at exempt
counts `0:sum(Zb)`.

## Details

Write \\r_0\\ for the ranks of the vector at exempt count zero, where
treated units sit at \\Y_i - c\\ and controls at \\Y_j\\. The null
exempts treated units in order of decreasing outcome, and every treated
unit shifts by the same \\c\\, so the treated units' order under \\r_0\\
is their outcome order and the \\ii\\ exempted units are the \\ii\\
treated units with the largest \\r_0\\.

Two things follow. No remaining unit has an exempted unit below it, so
each remaining rank is its \\r_0\\ shifted up by \\ii\\. And the
exempted units occupy ranks \\1, \ldots, ii\\, all of them treated, so
they contribute `sum(score[1:ii])` however they are ordered among
themselves. Hence

\$\$\mathrm{min\\stat}(ii) = \sum\_{l \le ii} s_l + \sum\_{j \> ii}
s\_{ii + r_0(u_j)},\$\$

where \\s\\ is the score vector and \\u_1, \ldots, u\_{m_b}\\ are the
treated units ordered by \\r_0\\. The identity is checked against
`min_stat` on random designs, with ties and at exact breakpoints, in
`tests/testthat/test-min_stat-path.R`.
