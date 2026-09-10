# CMRSS 0.2.11

## New

- The Gurobi solvers now honour `options(CMRSS.gurobi.threads = n)`. Gurobi's
  default is to use every core on the machine for a single solve, which is
  right for a caller solving one problem at a time and wrong for one that is
  already running in parallel. Running the combined_stephenson_tests SRE
  simulation as a `foreach` loop over 14 workers, each solve took 14 threads,
  which put roughly 150 runnable threads on a 14-core machine: the load average
  reached 158 and the workers together drew 581 percent of CPU where 14 busy
  cores would be 1400 percent. They were queueing rather than solving.

  Only the caller knows whether it is already parallel, so the thread count is
  the caller's to set. Forked workers inherit R options, so one
  `options(CMRSS.gurobi.threads = 1)` before a parallel loop reaches every
  worker. The option is unset by default and nothing changes for anyone who
  does not set it.

  `Threads` is a resource setting rather than a modelling one, so capping it
  does not change the answer. `tests/testthat/test-gurobi-threads.R` checks
  that on a four-block problem, comparing the bounds from a capped solve
  against an uncapped one from the same seed.

# CMRSS 0.2.10

## New exported functions

- Six building blocks that were internal are now part of the public
  interface: `sort_treat()`, `rank_score()`, `min_stat()`, `null_dist()`,
  `null_dist_multiple()` and `comb_null_dist_cre()`. The package already
  exported the functions that answer a whole question, such as a p-value or a
  set of confidence bounds. Anyone assembling a procedure of their own needed
  the pieces those are built from, and the only route was `CMRSS:::`. The
  combined_stephenson_tests paper repository had been carrying a 1,731-line
  copy of this package's code for exactly that reason. Each of the six now has
  a runnable example; none of their behaviour changed.

## Bug fixes

- `comb_null_dist_cre()` took the number of draws from its `nperm` argument
  while filling its tail-probability matrix from a null distribution that was
  `ncol(Z.perm)` wide. When the two disagreed R recycled the values to fit,
  and the result was a null distribution of the wrong length built from
  repeated entries. The count now comes from the null distribution in hand,
  which is what `null_dist_multiple()` and `com_conf_quant_larger_cre()`
  already did.

  This was reachable from outside the package through `comb_p_val_cre()`,
  which passes `nperm` straight down. On a 20-unit example with a fixed
  200-column `Z.perm`, `nperm = 200` gave `p = 0.8` and `nperm = 1000` gave
  `p = 0.000000`: the same data and the same permutations, and a rejection at
  any level, from an argument a caller could reasonably think was redundant
  once `Z.perm` was supplied. Callers who passed a matching `nperm`, including
  everything inside this package, are unaffected;
  `com_conf_quant_larger_cre()` has set `nperm <- ncol(Z.perm)` on entry all
  along.

## Tests

- `tests/testthat/test-exports.R` pins the public interface and checks that
  the six newly exported functions compute what they did when internal.
- `tests/testthat/test-nperm-zperm-consistency.R` checks that a supplied
  permutation matrix decides the number of draws, for both
  `comb_null_dist_cre()` and `comb_p_val_cre()`.

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
