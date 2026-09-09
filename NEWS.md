# CMRSS 0.3.0

## Confidence limits are now exact

- `com_block_conf_quant_larger()` and `com_block_conf_quant_larger_trt()`
  locate each limit by binary search over the thresholds at which the
  combined statistic can change, rather than by `uniroot()` followed by a
  walk in steps of `tol`. The statistic is a step function of the threshold
  `c` whose jumps occur only at the within-block differences `Y_i - Y_j`
  between a treated unit `i` and a control unit `j`, so the limit is one of
  those differences and can be found exactly.

- **This moves published numbers.** Limits shift by less than `tol` and
  always downward, because the old search could stop above the true
  infimum and so exclude thresholds the test does not reject. On the test
  design the four finite limits moved by -0.016, -0.038, -0.039, and
  -0.044 at `tol = 0.05`. Intervals are therefore slightly wider, and
  `tol` no longer affects the answer.

- New internal helpers in `R/breakpoints.R`: `min_stat_breakpoints()`,
  `interval_probes()`, `binary_search_first()`, `search_limit()`, and
  `memoize_on_threshold()`.

- The search evaluates the statistic strictly between breakpoints, never at
  one. At a breakpoint `Y_i - c` equals `Y_j` and `ties.method = "first"`
  resolves the tie by input order, so the value there follows unit ordering
  rather than the data. Keeping evaluations inside the intervals removes
  that dependence and leaves the search unaffected by PLAN item 2A.

- Coefficient matrices are cached on the threshold, since they do not
  depend on the quantile index while the optimization that consumes them
  does.

# CMRSS 0.2.9

## PLAN item 1B settled: the wider column range was redundant, not wrong

- `comb_matrix_block_stratum()` enumerated exempt counts `0:n_b` where the
  index counts treated units whose effect exceeds the threshold and so
  cannot exceed `m_b`. `min_stat()` exempts `min(m, n - k)` units, so every
  column past `m_b` repeated the `m_b` column's value while carrying a
  larger index into the stratum solver's budget constraint. Such a column
  uses more of that budget for the same objective, so no optimal solution
  used one.

- **No number changes.** The range is now `0:m_b`. `comb.method = 2`
  p-values are identical. The tests rebuild the wider matrix and confirm
  the solver returns the same objective at every budget from 0 to
  `sum(m_b)`.

- `comb_matrix_block_stratum()` also uses `min_stat_path()` now, and
  `max_comb_matrix_block_stratum()` forwards the caller's precomputed
  scores and block summary instead of passing `NULL` and making them be
  rebuilt at every threshold.

# CMRSS 0.2.8

## One ranking per block instead of m_b + 1

- `comb_matrix_block()` obtains the whole sequence of exempt counts from a
  single call to `rank()` through the new internal `min_stat_path()`, in
  place of `m_b + 1` separate `min_stat()` calls that each re-ran
  `sort_treat()` and `rank()` over the same block.

- **No number changes.** The identity is exact, not an approximation: the
  exempted units are the treated units with the largest ranks, so every
  remaining rank is its original shifted up by the exempt count, and the
  exempted units occupy the lowest ranks and are all treated, so they
  contribute a fixed partial sum of the scores. Checked against
  `min_stat()` on 1,200 random designs with ties and at exact breakpoints,
  agreeing to 1e-10, and `comb_matrix_block()` is pinned against stored
  values.

# CMRSS 0.2.7

## Documentation

- Clarified the treated-only convention of `pval_comb_block()` in the
  function's docstring and example: `k` is in `1..sum(Z)` (number of treated
  units), and the function tests the treated-only hypothesis
  `H_{k,c}^treat`. The example now uses `k <- floor(0.9 * sum(Z))`.

## Tests

- Re-enabled 10 tests previously skipped pending the k-convention resolution
  (`tests/testthat/test-CMRSS_SRE.R`, `test-pval_scre.R`, `test-pval-cre.R`,
  `test-solvers.R`). Each `k` formula was updated from the all-units
  convention to the treated-only convention.
- The two cross-validation tests against `RIQITE::pval_quantile`
  (`test-pval-cre.R`, `test-pval_scre.R`) now pass a shared `Z.perm`
  permutation matrix to both packages and apply the translation
  `k_R = k_C + (n - m)` between RIQITE's all-units `k` and CMRSS's
  treated-only `k`. Under this setup the two packages' p-values agree
  exactly.
- The four gurobi-gated solver-equivalence tests in `test-solvers.R` had the
  same latent k-convention bug; their `k` formulas are now treated-only.

# CMRSS 0.2.6

## Bug fixes

- `pval_comb_block()` now validates that `k` lies in `1..sum(Z)` and errors
  with a clear message otherwise. Previously, `k > sum(Z)` made the LP
  infeasible and the function silently returned `p.value = 0` with
  `test.stat = Inf` (a false rejection). The function tests the treated-only
  hypothesis `H_{k,c}^treat`, so `k` must not exceed the number of treated
  units.

# CMRSS 0.2.5

## Performance

- Avoid large `n x nperm` permutation matrices by default:
  - CRE null generation (`null_dist()`, `null_dist_multiple()`) now streams permutations in chunks when `Z.perm` is not supplied (new `chunk_size` argument).
  - SRE null generation (`com_null_dist_block()`, `com_null_dist_block_stratum()`) now generates permutations in chunks when `Z.perm` is not supplied (new `chunk_size` argument).
- Speed up stratified runs with many blocks:
  - `summary_block()` now builds `units.block` via `split()` (reduces `O(n·B)` scanning to `O(n)`).
  - Block-randomized chunk generation uses fast paths for common small-block cases (e.g., `mb == 1`, `mb == nb - 1`) and an internal cache for small `combn(nb, mb)` patterns.
- Speed up `comb.method = 2` null distribution:
  - `com_null_dist_block_stratum()` now computes per-block statistics with matrix operations instead of per-permutation nested loops.
- Reduce allocation/copying overhead in solver setup:
  - Solver constraint triplets are now assembled via list accumulation rather than repeated `c()` growth.
- Minor speed improvement in generalized CI code:
  - `ci_lower_quantile_generalize()` no longer grows results with repeated `rbind()` inside a loop.

## Bug Fixes

- `rank_score()` consistently returns a numeric vector (avoids 1-column matrix output when `scale = TRUE`).

## Testing

- Added unit tests for internal block-permutation chunk generation and for equivalence of the vectorized `com_null_dist_block_stratum()` implementation to a naive reference computation.
