# Changelog

## CMRSS 0.2.16

### Documentation

- `README.md` installed the package from `jwbowers/CMRSS` and showed
  that repository’s check badge. GitHub now redirects `jwbowers/CMRSS`
  to the archived `bowers-illinois-edu/CMRSS_archive`, so the
  instructions installed an old version. The README now installs from
  `davidk91919/CMRSS`, shows that repository’s badge, and links the
  pkgdown site. `DESCRIPTION` points `URL` and `BugReports` at
  `davidk91919/CMRSS`, where the issues are filed.

## CMRSS 0.2.15

### Documentation

- New vignette, “Could anyone have been harmed?”
  ([`vignette("harm", package = "CMRSS")`](https://bowers-illinois-edu.github.io/CMRSS/articles/harm.md)).
  It shows how analyzing minus the outcome turns the package’s lower
  bounds into a lower count of harmed units, how the count of helped
  units gives an upper count of units that could have been harmed, and
  why the two analyses run at `alpha / 2` when reported together. It
  works through a simulated experiment in which 40 of 200 units are
  harmed while the average effect is positive, and through
  `electric_teachers`. `knitr` and `rmarkdown` join Suggests.

### Other changes

- The HiGHS solvers no longer print “LP matrix packed vector contains …
  values … ignored” on every solve. Stephenson scores with a large s
  give constraint-matrix entries at or below 1e-9, which HiGHS ignores
  and warns about. The package now removes them before calling HiGHS, so
  HiGHS solves the same problem as before. On `electric_teachers`,
  within sites, with s = 2, 11, 59, every bound from
  [`cmrss()`](https://bowers-illinois-edu.github.io/CMRSS/reference/cmrss.md)
  is identical before and after.

## CMRSS 0.2.14

### New

- [`cmrss()`](https://bowers-illinois-edu.github.io/CMRSS/reference/cmrss.md)
  reads an experiment from a formula and a data frame: `gain ~ TxAny`
  for a completely randomized experiment, `gain ~ TxAny | Site` for a
  block-randomized one. It returns lower confidence bounds for every
  sorted effect in the chosen set (`"treat"`, `"control"` or `"all"`)
  and, given `quantile`, a p-value for the hypothesis that the effect at
  that proportion of the set is at most `c`. It converts the proportion
  into the `k` each underlying function expects, which differs between
  [`comb_p_val_cre()`](https://bowers-illinois-edu.github.io/CMRSS/reference/comb_p_val_cre.md)
  (all units) and
  [`pval_comb_block()`](https://bowers-illinois-edu.github.io/CMRSS/reference/pval_comb_block.md)
  (treated units), and reports the `k` used. Writing `-gain` on the left
  bounds effects on minus the outcome, which is how to ask whether
  anyone could have been harmed.
  [`print()`](https://rdrr.io/r/base/print.html) reports how many units
  have bounds above `c`.
- The rank statistics are polynomial scores, (r/(n + 1))^(s - 1), with
  parameters `s` (the paper’s zeta); `scores = "stephenson"` uses
  Stephenson scores choose(r - 1, s - 1) instead. The default is three
  values, 2, s_max and round(sqrt(2 s_max)), where s_max = 4 m / q_min
  (at most n / 2) and q_min is the fewest treated units whose holding
  the top ranks would give a p-value at or below `alpha`. The section
  “Choosing s” of
  [`?cmrss`](https://bowers-illinois-edu.github.io/CMRSS/reference/cmrss.md)
  gives the reasoning and a simulation comparing grids. Polynomial
  scores are never 0, so they work in blocks of any size; with
  Stephenson scores, an s above a block’s size is lowered to that size
  so the block is not dropped.
- [`comb_p_val_cre()`](https://bowers-illinois-edu.github.io/CMRSS/reference/comb_p_val_cre.md),
  [`com_conf_quant_larger_cre()`](https://bowers-illinois-edu.github.io/CMRSS/reference/com_conf_quant_larger_cre.md),
  [`pval_comb_block()`](https://bowers-illinois-edu.github.io/CMRSS/reference/pval_comb_block.md),
  [`com_block_conf_quant_larger()`](https://bowers-illinois-edu.github.io/CMRSS/reference/com_block_conf_quant_larger.md)
  and
  [`cmrss()`](https://bowers-illinois-edu.github.io/CMRSS/reference/cmrss.md)
  now warn, with class `"cmrss_ties_warning"`, when the outcome takes
  two values or when one value is shared by at least two units and more
  than 5 percent of them. Tied outcomes are ranked by row order, so
  reordering rows can change results (issue
  [\#5](https://github.com/davidk91919/CMRSS/issues/5)).
  `electric_teachers$gain` triggers the warning: 29 of 233 teachers
  share one value. No returned number changes.

## CMRSS 0.2.13

### Bug fixes

- The two stratified examples in `README.md`, the Quick Start and
  Example 1, stopped with an error. They passed `k = floor(0.9 * N)`,
  counted over all N units, to
  [`pval_comb_block()`](https://bowers-illinois-edu.github.io/CMRSS/reference/pval_comb_block.md),
  which since 0.2.7 counts `k` over the treated units only and so
  requires `1 <= k <= sum(Z)`. Both now use `k = floor(0.9 * sum(Z))`.
  The “Comparing Solvers” example used `k` and `c` without defining
  them; the Quick Start now defines both, and the comparison sets the
  same seed before each solver so the two null distributions match.
- [`pval_comb_block()`](https://bowers-illinois-edu.github.io/CMRSS/reference/pval_comb_block.md)’s
  error for an out-of-range `k` now says that `k` counts treated units,
  gives the range `1..sum(Z)`, and points to
  `com_block_conf_quant_larger(set = "all")` for quantiles of all units.
  It used to cite a source line, `R/CMRSS_SRE.R:1034`, that a user of
  the installed package cannot see.
- [`pval_comb_block()`](https://bowers-illinois-edu.github.io/CMRSS/reference/pval_comb_block.md)
  now refuses a `k` that is not a whole number. It used to return a
  p-value for, say, `k = 2.5`.

`tests/testthat/test-readme-examples.R` runs every R chunk of the README
from a source checkout, so a stale example now fails `devtools::test()`.
No returned number changes.

## CMRSS 0.2.12

### Bug fixes

- [`com_block_conf_quant_larger()`](https://bowers-illinois-edu.github.io/CMRSS/reference/com_block_conf_quant_larger.md)
  now relabels a caller-supplied `Z.perm` for the control bounds. Those
  bounds come from the relabeled experiment, in which `Z` becomes
  `1 - Z` and `Y` becomes `-Y`, and each simulated assignment has to be
  relabeled the same way. The function passed `Z.perm` through unchanged
  for `set = "control"` and for the control half of `set = "all"`, so
  every simulated assignment scored the wrong units; on the
  `electric_teachers` data every control bound came back `-Inf`.
  [`com_conf_quant_larger_cre()`](https://bowers-illinois-edu.github.io/CMRSS/reference/com_conf_quant_larger_cre.md)
  already relabeled `Z.perm`.
- The stratified interval function now counts the simulated assignments
  it actually has when it picks the critical value. It took position
  `floor(null.max * alpha) + 1` of the sorted null distribution using
  the argument `null.max`, even when the caller supplied a `Z.perm` or
  `stat.null` of a different length. With a 2000-column `Z.perm`, the
  default `null.max = 10^4` and `alpha = 0.10`, it took position 1001 of
  2000 and so tested at level 0.50. `null.max` is now set to
  `length(stat.null)`.
- `com_block_conf_quant_larger(set = "all")` now stops when given
  `stat.null`. The two halves of `set = "all"` need different null
  distributions whenever a stratum is not half treated, so one vector
  cannot serve both.

Results change only for callers who pass `Z.perm` or `stat.null`. When
neither is passed, `length(stat.null)` equals `null.max` and the bounds
are unchanged. None of the stratified scripts in the
combined_stephenson_tests paper passes either argument. Tests are in
`tests/testthat/test-sre-ci-user-null.R`.

## CMRSS 0.2.11

### New

- The Gurobi solvers now honour `options(CMRSS.gurobi.threads = n)`.
  Gurobi’s default is to use every core on the machine for a single
  solve, which is right for a caller solving one problem at a time and
  wrong for one that is already running in parallel. Running the
  combined_stephenson_tests SRE simulation as a `foreach` loop over 14
  workers, each solve took 14 threads, which put roughly 150 runnable
  threads on a 14-core machine: the load average reached 158 and the
  workers together drew 581 percent of CPU where 14 busy cores would be
  1400 percent. They were queueing rather than solving.

  Only the caller knows whether it is already parallel, so the thread
  count is the caller’s to set. Forked workers inherit R options, so one
  `options(CMRSS.gurobi.threads = 1)` before a parallel loop reaches
  every worker. The option is unset by default and nothing changes for
  anyone who does not set it.

  `Threads` is a resource setting rather than a modelling one, so
  capping it does not change the answer.
  `tests/testthat/test-gurobi-threads.R` checks that on a four-block
  problem, comparing the bounds from a capped solve against an uncapped
  one from the same seed.

## CMRSS 0.2.10

### New exported functions

- Six building blocks that were internal are now part of the public
  interface:
  [`sort_treat()`](https://bowers-illinois-edu.github.io/CMRSS/reference/sort_treat.md),
  [`rank_score()`](https://bowers-illinois-edu.github.io/CMRSS/reference/rank_score.md),
  [`min_stat()`](https://bowers-illinois-edu.github.io/CMRSS/reference/min_stat.md),
  [`null_dist()`](https://bowers-illinois-edu.github.io/CMRSS/reference/null_dist.md),
  [`null_dist_multiple()`](https://bowers-illinois-edu.github.io/CMRSS/reference/null_dist_multiple.md)
  and
  [`comb_null_dist_cre()`](https://bowers-illinois-edu.github.io/CMRSS/reference/comb_null_dist_cre.md).
  The package already exported the functions that answer a whole
  question, such as a p-value or a set of confidence bounds. Anyone
  assembling a procedure of their own needed the pieces those are built
  from, and the only route was `CMRSS:::`. The combined_stephenson_tests
  paper repository had been carrying a 1,731-line copy of this package’s
  code for exactly that reason. Each of the six now has a runnable
  example; none of their behaviour changed.

### Bug fixes

- [`comb_null_dist_cre()`](https://bowers-illinois-edu.github.io/CMRSS/reference/comb_null_dist_cre.md)
  took the number of draws from its `nperm` argument while filling its
  tail-probability matrix from a null distribution that was
  `ncol(Z.perm)` wide. When the two disagreed R recycled the values to
  fit, and the result was a null distribution of the wrong length built
  from repeated entries. The count now comes from the null distribution
  in hand, which is what
  [`null_dist_multiple()`](https://bowers-illinois-edu.github.io/CMRSS/reference/null_dist_multiple.md)
  and
  [`com_conf_quant_larger_cre()`](https://bowers-illinois-edu.github.io/CMRSS/reference/com_conf_quant_larger_cre.md)
  already did.

  This was reachable from outside the package through
  [`comb_p_val_cre()`](https://bowers-illinois-edu.github.io/CMRSS/reference/comb_p_val_cre.md),
  which passes `nperm` straight down. On a 20-unit example with a fixed
  200-column `Z.perm`, `nperm = 200` gave `p = 0.8` and `nperm = 1000`
  gave `p = 0.000000`: the same data and the same permutations, and a
  rejection at any level, from an argument a caller could reasonably
  think was redundant once `Z.perm` was supplied. Callers who passed a
  matching `nperm`, including everything inside this package, are
  unaffected;
  [`com_conf_quant_larger_cre()`](https://bowers-illinois-edu.github.io/CMRSS/reference/com_conf_quant_larger_cre.md)
  has set `nperm <- ncol(Z.perm)` on entry all along.

### Tests

- `tests/testthat/test-exports.R` pins the public interface and checks
  that the six newly exported functions compute what they did when
  internal.
- `tests/testthat/test-nperm-zperm-consistency.R` checks that a supplied
  permutation matrix decides the number of draws, for both
  [`comb_null_dist_cre()`](https://bowers-illinois-edu.github.io/CMRSS/reference/comb_null_dist_cre.md)
  and
  [`comb_p_val_cre()`](https://bowers-illinois-edu.github.io/CMRSS/reference/comb_p_val_cre.md).

## CMRSS 0.2.9

### PLAN item 1B settled: the wider column range was redundant, not wrong

- [`comb_matrix_block_stratum()`](https://bowers-illinois-edu.github.io/CMRSS/reference/comb_matrix_block_stratum.md)
  enumerated exempt counts `0:n_b` where the index counts treated units
  whose effect exceeds the threshold and so cannot exceed `m_b`.
  [`min_stat()`](https://bowers-illinois-edu.github.io/CMRSS/reference/min_stat.md)
  exempts `min(m, n - k)` units, so every column past `m_b` repeated the
  `m_b` column’s value while carrying a larger index into the stratum
  solver’s budget constraint. Such a column uses more of that budget for
  the same objective, so no optimal solution used one.

- **No number changes.** The range is now `0:m_b`. `comb.method = 2`
  p-values are identical. The tests rebuild the wider matrix and confirm
  the solver returns the same objective at every budget from 0 to
  `sum(m_b)`.

- [`comb_matrix_block_stratum()`](https://bowers-illinois-edu.github.io/CMRSS/reference/comb_matrix_block_stratum.md)
  also uses
  [`min_stat_path()`](https://bowers-illinois-edu.github.io/CMRSS/reference/min_stat_path.md)
  now, and
  [`max_comb_matrix_block_stratum()`](https://bowers-illinois-edu.github.io/CMRSS/reference/max_comb_matrix_block_stratum.md)
  forwards the caller’s precomputed scores and block summary instead of
  passing `NULL` and making them be rebuilt at every threshold.

## CMRSS 0.2.8

### One ranking per block instead of m_b + 1

- [`comb_matrix_block()`](https://bowers-illinois-edu.github.io/CMRSS/reference/comb_matrix_block.md)
  obtains the whole sequence of exempt counts from a single call to
  [`rank()`](https://rdrr.io/r/base/rank.html) through the new internal
  [`min_stat_path()`](https://bowers-illinois-edu.github.io/CMRSS/reference/min_stat_path.md),
  in place of `m_b + 1` separate
  [`min_stat()`](https://bowers-illinois-edu.github.io/CMRSS/reference/min_stat.md)
  calls that each re-ran
  [`sort_treat()`](https://bowers-illinois-edu.github.io/CMRSS/reference/sort_treat.md)
  and [`rank()`](https://rdrr.io/r/base/rank.html) over the same block.

- **No number changes.** The identity is exact, not an approximation:
  the exempted units are the treated units with the largest ranks, so
  every remaining rank is its original shifted up by the exempt count,
  and the exempted units occupy the lowest ranks and are all treated, so
  they contribute a fixed partial sum of the scores. Checked against
  [`min_stat()`](https://bowers-illinois-edu.github.io/CMRSS/reference/min_stat.md)
  on 1,200 random designs with ties and at exact breakpoints, agreeing
  to 1e-10, and
  [`comb_matrix_block()`](https://bowers-illinois-edu.github.io/CMRSS/reference/comb_matrix_block.md)
  is pinned against stored values.

## CMRSS 0.2.7

### Documentation

- Clarified the treated-only convention of
  [`pval_comb_block()`](https://bowers-illinois-edu.github.io/CMRSS/reference/pval_comb_block.md)
  in the function’s docstring and example: `k` is in `1..sum(Z)` (number
  of treated units), and the function tests the treated-only hypothesis
  `H_{k,c}^treat`. The example now uses `k <- floor(0.9 * sum(Z))`.

### Tests

- Re-enabled 10 tests previously skipped pending the k-convention
  resolution (`tests/testthat/test-CMRSS_SRE.R`, `test-pval_scre.R`,
  `test-pval-cre.R`, `test-solvers.R`). Each `k` formula was updated
  from the all-units convention to the treated-only convention.
- The two cross-validation tests against
  [`RIQITE::pval_quantile`](https://rdrr.io/pkg/RIQITE/man/pval_quantile.html)
  (`test-pval-cre.R`, `test-pval_scre.R`) now pass a shared `Z.perm`
  permutation matrix to both packages and apply the translation
  `k_R = k_C + (n - m)` between RIQITE’s all-units `k` and CMRSS’s
  treated-only `k`. Under this setup the two packages’ p-values agree
  exactly.
- The four gurobi-gated solver-equivalence tests in `test-solvers.R` had
  the same latent k-convention bug; their `k` formulas are now
  treated-only.

## CMRSS 0.2.6

### Bug fixes

- [`pval_comb_block()`](https://bowers-illinois-edu.github.io/CMRSS/reference/pval_comb_block.md)
  now validates that `k` lies in `1..sum(Z)` and errors with a clear
  message otherwise. Previously, `k > sum(Z)` made the LP infeasible and
  the function silently returned `p.value = 0` with `test.stat = Inf` (a
  false rejection). The function tests the treated-only hypothesis
  `H_{k,c}^treat`, so `k` must not exceed the number of treated units.

## CMRSS 0.2.5

### Performance

- Avoid large `n x nperm` permutation matrices by default:
  - CRE null generation
    ([`null_dist()`](https://bowers-illinois-edu.github.io/CMRSS/reference/null_dist.md),
    [`null_dist_multiple()`](https://bowers-illinois-edu.github.io/CMRSS/reference/null_dist_multiple.md))
    now streams permutations in chunks when `Z.perm` is not supplied
    (new `chunk_size` argument).
  - SRE null generation
    ([`com_null_dist_block()`](https://bowers-illinois-edu.github.io/CMRSS/reference/com_null_dist_block.md),
    [`com_null_dist_block_stratum()`](https://bowers-illinois-edu.github.io/CMRSS/reference/com_null_dist_block_stratum.md))
    now generates permutations in chunks when `Z.perm` is not supplied
    (new `chunk_size` argument).
- Speed up stratified runs with many blocks:
  - [`summary_block()`](https://bowers-illinois-edu.github.io/CMRSS/reference/summary_block.md)
    now builds `units.block` via
    [`split()`](https://rdrr.io/r/base/split.html) (reduces `O(n·B)`
    scanning to `O(n)`).
  - Block-randomized chunk generation uses fast paths for common
    small-block cases (e.g., `mb == 1`, `mb == nb - 1`) and an internal
    cache for small `combn(nb, mb)` patterns.
- Speed up `comb.method = 2` null distribution:
  - [`com_null_dist_block_stratum()`](https://bowers-illinois-edu.github.io/CMRSS/reference/com_null_dist_block_stratum.md)
    now computes per-block statistics with matrix operations instead of
    per-permutation nested loops.
- Reduce allocation/copying overhead in solver setup:
  - Solver constraint triplets are now assembled via list accumulation
    rather than repeated [`c()`](https://rdrr.io/r/base/c.html) growth.
- Minor speed improvement in generalized CI code:
  - [`ci_lower_quantile_generalize()`](https://bowers-illinois-edu.github.io/CMRSS/reference/ci_lower_quantile_generalize.md)
    no longer grows results with repeated
    [`rbind()`](https://rdrr.io/r/base/cbind.html) inside a loop.

### Bug Fixes

- [`rank_score()`](https://bowers-illinois-edu.github.io/CMRSS/reference/rank_score.md)
  consistently returns a numeric vector (avoids 1-column matrix output
  when `scale = TRUE`).

### Testing

- Added unit tests for internal block-permutation chunk generation and
  for equivalence of the vectorized
  [`com_null_dist_block_stratum()`](https://bowers-illinois-edu.github.io/CMRSS/reference/com_null_dist_block_stratum.md)
  implementation to a naive reference computation.
