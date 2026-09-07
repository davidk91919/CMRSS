# Changelog

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
