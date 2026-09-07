# CMRSS Audit Implementation Plan

A numbered checkbox plan grouped by priority. Each task lists:
investigation, tests-first, implementation, and a Jake review
checkpoint. Honors `~/repos/ai_workflow/CLAUDE_CODING.md`: tests before
code, pause after tests, pause after implementation, pause at any
unresolved design decision, bump DESCRIPTION patch version when
adding/exporting, comment WHY not WHAT.

Independence map (for parallelization):

- 1A and 1B both touch `R/CMRSS_SRE.R` and the same coefficient/LP path;
  do 1A first then 1B (or one branch each).
- 2A is independent of 1A/1B and can run in a separate worktree.
- 2B depends on 1A being resolved (the convention for `p` and `k` must
  be settled before flipping the default).
- 2C is fully independent (`R/comparison_methods.R`).
- 3A, 3B, 3E are pure comments and can be done together at the end.
- 3C and 3D touch user-facing defaults; bundle with 2B if Jake decides
  defaults all change at once.
- 3F is a follow-on test bundle that pins anything not already pinned by
  1A/1B/2B.

Suggested suites of parallel work (worktrees):

- Worktree A: 1A + 1B + 2B + 3F (the LP/`p`/`k` math)
- Worktree B: 2A (ties) + 2C (renaming Berger-Boos)
- Worktree C: 3A, 3B, 3C, 3D, 3E (docs/cosmetics) after A and B land

Version bumping rule used below: any change that adds, removes, or
renames an exported symbol bumps DESCRIPTION patch (currently 0.2.5).
Internal-only fixes still bump patch when they alter user-visible
numerical results (1A, 1B, 2A, 2B).

------------------------------------------------------------------------

## PRIORITY 1 – Possible bugs (investigate; pin with tests; fix if confirmed)

### \[x\] 1A. Reconcile the `p` convention in `pval_comb_block` vs `com_block_conf_quant_larger_trt`

**Closed 2026-05-19 (treated-only convention; tests re-enabled in CMRSS
0.2.7).** Resolution: `R/CMRSS_SRE.R:1034` (`p <- m - k`) is correct.
`pval_comb_block` tests the treated-only hypothesis `H_{k,c}^treat` with
`k in 1..sum(Z)`. The docstring and example have been updated (commit
forthcoming; see NEWS.md 0.2.7). The 10 tests previously skipped pending
this resolution have been re-enabled with `k` formulas converted from
the all-units convention to the treated-only convention. Four additional
gurobi-gated sites in `test-solvers.R` had the same latent bug and were
fixed in the same pass. The two cross-validation tests against
[`RIQITE::pval_quantile`](https://rdrr.io/pkg/RIQITE/man/pval_quantile.html)
now pass a shared `Z.perm` to both packages and apply the LP-equivalence
translation `k_R = k_C + (n - m)`; under that setup the p-values agree
exactly. The input-validation guard added at `R/CMRSS_SRE.R:1040` in
0.2.6 stays in place. `devtools::check()` is clean (0 errors, 0
warnings, 2 pre-existing repo-hygiene NOTEs unrelated to this item).

Still open from the larger scenario-(a) audit:

- `com_block_conf_quant_larger` `set = "all"` and `set = "control"`
  paths (`R/CMRSS_SRE.R:1430-1499`): **Closed 2026-05-19, no functional
  bug.** Verified against paper eq:H_kc_t / eq:H_kc_c (main.tex lines
  459-466, 479, 502, 940, 948-950). The swap `Z <- 1 - Z, Y <- -Y`
  preserves individual effects via the relabeling
  `tilde_tau_treat(k) = tau_control(k)` (paper line 479), so the
  swap-space lower CIs are directly lower CIs on `tau_control(k)` in the
  original problem; no sign flip or index reversal is needed.
  `set = "all"` correctly applies the Bonferroni split (alpha/2 per
  branch) per paper line 502. New tests in
  `tests/testthat/test-com_block_conf_quant_larger-sets.R` pin the
  wrapper’s structural behavior at `tolerance = 1e-6` via paired
  `set.seed` calls. `com_block_conf_quant_larger_trt`’s function
  documentation expanded to record the all-units convention and the
  LP-equivalence with `pval_comb_block`.
- `com_block_conf_quant_larger_trt` (`R/CMRSS_SRE.R:1171, 1234`) still
  uses `p <- n - k` (all-units). David did not touch these in
  `762c4d08`. Now correctly *documented* as the all-units convention;
  remaining question is whether to refactor `_trt` to take treated-only
  `k` internally (functionally equivalent, removes the asymmetry with
  `pval_comb_block`, but risks breaking any downstream caller of
  `_trt`). Deferred.
- Paper experiment scripts in
  `/Users/jwbowers/repos/combined_stephenson_tests/code/`: with the
  wrapper now confirmed correct, paper-side numbers are unaffected. Low
  priority; deferred.

Item 2B (default `comb.method` flip) is now unblocked.

Status update from initial reproducer:
`pval_comb_block(Z, Y, k=19, c=0, ...)` with the existing test inputs
(s=3, n=8, m=4 per stratum, total m=12) yields `p = m - k = -7`, the LP
returns “Infeasible,” the test statistic is `Inf`, and
`mean(stat.null >= Inf) = 0`. So the function silently returns
`p.value = 0` whenever k \> m_total. This is a real, consequential bug.

1.  Investigate

Re-derive the LP coverage constraint from the paper at
`/Users/jwbowers/repos/combined_stephenson_tests/main.tex`
(eq:constraint_ks). Determine whether `p` should be `m - k`
(treated-only) or `n - k` (over all units), as a function of which `k`
regime is being tested (effect-quantile-among-treated vs all-quantile).

Read `/Users/jwbowers/repos/CMRSS_jake/R/CMRSS_SRE.R` lines 1015-1085
(`pval_comb_block`, line 1034 sets `p <- m - k`) and lines 1102-1289
(`com_block_conf_quant_larger_trt`, line 1171 and line 1234 set
`p <- n - k`).

Trace the wrapper `com_block_conf_quant_larger` (lines 1416-1485) to see
that `set = "treat"` calls `_trt` directly (so `_trt` operates on the
original `Z`), but `set = "control"` and `set = "all"` call it after
swapping `Z <- 1 - Z`. After the swap, `m` and `n - m` switch roles;
this is where the divergence between `m - k` and `n - k` matters most.

Read `/Users/jwbowers/repos/CMRSS_jake/tests/testthat/test-CMRSS_SRE.R`
(esp. tests using `k = floor(0.8 * N)`) to confirm whether they
currently exercise both code paths and whether negative `p` is ever
silently produced.

Read `solve_optimization` (`R/solvers.R` lines 422-447) and
`HiGHS_sol_com` (lines 104-252): find what happens when `p` is negative
or larger than `sum(mb)`. The coverage row is `sum(k_x * x) <= p`; with
negative `p` the LP is infeasible. Confirm whether HiGHS returns status
!= 7 (warning at line 248) or silently returns `Inf`/`NA`.

Build a tiny brute-force reproducer outside the package: B = 2 blocks,
`nb = c(4, 3)`, `mb = c(2, 1)`, fixed Y, two methods, all `2^N`
xi-allocations consistent with each candidate `p`, take the min
standardized combined statistic. Compare against `solve_optimization`
output for both `p = m - k` and `p = n - k` formulas.

Files: `R/CMRSS_SRE.R:1034`, `R/CMRSS_SRE.R:1171`, `R/CMRSS_SRE.R:1234`,
`R/solvers.R:104-252`, `R/solvers.R:422-447`,
`tests/testthat/test-CMRSS_SRE.R`.

2.  Tests-first (write before any code change; commit; pause)

Add `tests/testthat/test-pval-comb-block-p-convention.R` containing:

- A brute-force enumeration helper for tiny SRE problems that returns
  the true min combined statistic given a candidate `p`.
- One
  `test_that("pval_comb_block p convention matches paper eq:constraint_ks ...")`
  that calls `pval_comb_block` with k corresponding to the
  treated-quantile regime and compares to brute force. This pins what
  `p` should be.
- One
  `test_that("com_block_conf_quant_larger_trt p convention is consistent with pval_comb_block")`
  that runs both functions on the same tiny problem and asserts the LP
  solves are equal (after accounting for any documented sign /
  direction).
- One `test_that("infeasible p is detected, not silently solved")` that
  constructs `m`, `k` so that `m - k < 0` and expects either an explicit
  error or a documented sentinel value (Jake decides).

Run tests, confirm at least one fails as expected. Pause.

CHECKPOINT 1A.t (after tests): Jake reviews tests, confirms which
convention is correct, and decides the behavior for infeasible `p`.

3.  Implementation

Harmonize the two functions to the agreed convention. If the agreed
convention is `p = n - k` everywhere, change `R/CMRSS_SRE.R:1034`. If it
is `p = m - k`, change `R/CMRSS_SRE.R:1171` and `R/CMRSS_SRE.R:1234`.
Add a WHY comment that cites the paper equation.

Add an explicit feasibility guard near the top of `pval_comb_block` and
`com_block_conf_quant_larger_trt` that errors (or returns a documented
sentinel) when `p < 0` or `p > sum(mb)`. WHY comment: silent LP
infeasibility produced misleading output prior to this guard.

Bump DESCRIPTION patch (0.2.5 -\> 0.2.6) since user-visible numerical
output may change.

Run `devtools::document()` and `devtools::check()`.

CHECKPOINT 1A.i (after implementation): Jake reviews diff and
`devtools::check()` output before merge.

------------------------------------------------------------------------

### \[ \] 1B. Verify column range `0:nb[i]` in `comb_matrix_block_stratum` against paper eq:comb_per_stratum

1.  Investigate

Compare `comb_matrix_block` (`R/CMRSS_SRE.R:487-535`, line 520 uses
`0:mb[i]`) to `comb_matrix_block_stratum` (`R/CMRSS_SRE.R:661-719`, line
704 uses `0:nb[i]`). They use different ranges for what looks like the
same enumeration over “number of treated units in this block whose
effect exceeds c.”

Read paper eq:comb_per_stratum (search
`combined_stephenson_tests/main.tex` for `comb_per_stratum`) and check
whether the index runs over treated count (`0..m_b`) or block size
(`0..n_b`).

Trace how `coeflist` from `comb_matrix_block_stratum` flows into
`solve_stratum_optimization` (`R/solvers.R:662-686` and
`Gurobi_sol_stratum_com:523-560`, `HiGHS_sol_stratum_com:577-642`). The
LP variable `x_{s,j}` is read off `coeflist[[i]][1, ]` (the row of
indices). If indices range up to `nb[i]` but only `0..mb[i]` are
meaningful, the extras either get pruned by other constraints or
contribute spurious feasible solutions.

Files: `R/CMRSS_SRE.R:487-535`, `R/CMRSS_SRE.R:661-719`,
`R/solvers.R:523-686`.

2.  Tests-first

Add `tests/testthat/test-comb-matrix-block-stratum-range.R`:

- On the tiny SRE example from 1A, brute-force the correct stratum-wise
  standardized objective values and the set of feasible `j` indices.
- `test_that("comb_matrix_block_stratum index range matches paper eq:comb_per_stratum")`
  that compares the returned `Ti[1, ]` row against the brute-force
  feasible index set.
- `test_that("comb.method = 1 and comb.method = 2 give same min value when correctly aligned")`
  on the same tiny example, where both code paths are valid; this
  catches a mismatch where one over-enumerates.

Run, confirm failure. Pause.

CHECKPOINT 1B.t: Jake reviews and decides whether to tighten to
`0:mb[i]` or document why `0:nb[i]` is intentional (e.g., to allow
xi-augmented control units).

3.  Implementation

Either change line 704 to `0 : mb[i]` and the loop on line 709
accordingly, OR add a WHY comment explaining the broader range. If
tightening, also confirm `solve_stratum_optimization` is unaffected (it
already reads `ncol(coeflist[[i]])`, so it adapts).

Bump DESCRIPTION patch if numerical outputs change.

CHECKPOINT 1B.i: Jake reviews diff and `devtools::check()`.

------------------------------------------------------------------------

## PRIORITY 2 – Paper-code consistency

### \[ \] 2A. Tie handling in `rank()` calls

Note: a quick scan suggests the existing code uses
`ties.method = "first"` (not `"average"` as the audit memo states):
`R/CMRSS_CRE.R:13, 278`; `R/CMRSS_SRE.R:231, 799`. The validity question
is the same: `"first"` is deterministic and not the random-tie-breaking
rule that finite-sample distribution-free validity sometimes requires.
Verify the actual `ties.method` argument before acting.

1.  Investigate

Confirm by reading `main.tex` (search for “tie”, “ties”, “ties.method”,
“random”) whether the paper’s validity proof assumes random
tie-breaking, or whether the rank-sum statistic is invariant to the tie
rule under the paper’s hypotheses.

Check the `electric_teachers` outcome (`R/data.R`, `data/`) for actual
ties; quantify how often ties occur in the worked examples.

Decide whether the right move is: (i) change default to `"random"` (Jake
then needs to set seeds for reproducibility), (ii) keep the current
default and document the assumption, or (iii) expose `ties.method` as a
function argument.

2.  Tests-first

Add `tests/testthat/test-tie-handling.R`:

- Construct a CRE example with heavy ties (binary or rounded outcomes),
  fixed Z, run `comb_p_val_cre` many times with different random seeds.
  Assert the desired invariance (or, if `"random"`, assert the empirical
  type I error under a true null is approximately alpha across seeds).
- Same for an SRE example via `pval_comb_block`.

Run, document current behavior. Pause.

CHECKPOINT 2A.t: Jake decides between options (i)/(ii)/(iii).

3.  Implementation (one of)

\(i\) Change defaults to `"random"` in the relevant
[`rank()`](https://rdrr.io/r/base/rank.html) calls. Document that
callers needing reproducibility should
[`set.seed()`](https://rdrr.io/r/base/Random.html). WHY comment cites
finite-sample validity.

\(ii\) Add a roxygen `@details` block to `pval_cre`, `min_stat`,
`comb_p_val_cre`, `pval_comb_block` explaining why the current default
is acceptable.

\(iii\) Add a `ties.method` argument plumbed through to all relevant
functions, default to current, advertise `"random"` for ties-heavy data.
Bump DESCRIPTION patch (new exported argument).

CHECKPOINT 2A.i: review.

------------------------------------------------------------------------

### \[ \] 2B. Flip default `comb.method` from 1 to 2

1.  Investigate

Confirm by reading `combined_stephenson_tests/code/` example scripts
that they use `comb.method = 2`. Confirm Section 4 of `main.tex`
recommends method 2 as the Theorem 5 setting.

Find all callers of `pval_comb_block` and `com_block_conf_quant_larger`:
Grep across `/Users/jwbowers/repos/CMRSS_jake` and
`/Users/jwbowers/repos/combined_stephenson_tests`. Note any caller that
relies on the implicit default.

Read `R/CMRSS_SRE.R:1023` (`comb.method = 1` default in
`pval_comb_block`) and `R/CMRSS_SRE.R:1422` (`comb.method = 1` default
in `com_block_conf_quant_larger`).

2.  Tests-first

Add `tests/testthat/test-comb-method-default.R`:

- `test_that("default comb.method is 2")` calls `pval_comb_block` and
  `com_block_conf_quant_larger` without specifying `comb.method` and
  asserts behavior matches an explicit `comb.method = 2` call (e.g.,
  same numerical output on a fixed seed).
- `test_that("comb.method = 2 is at least as powerful as 1 on a Theorem 5 case")`
  runs both on a synthetic example with a tail effect and asserts the
  rejection rate (over many seeds) for method 2 \>= method 1, or asserts
  the lower confidence limit for method 2 \>= method 1 (more
  informative). This pins the substantive justification, not just the
  default.

Run, observe failure. Pause.

CHECKPOINT 2B.t: Jake confirms the flip and the power claim under
Theorem 5.

3.  Implementation

Change default `comb.method = 1` to `comb.method = 2` at
`R/CMRSS_SRE.R:1023` and `R/CMRSS_SRE.R:1422`.

Update roxygen docs: in the `@param comb.method` blocks
(`R/CMRSS_SRE.R:968-973` and `R/CMRSS_SRE.R:1326-1330`), mark “2” as
default and note that method 1 is retained for reproducing earlier
results.

Update existing tests that call without `comb.method` to either pass
`comb.method = 1` explicitly (if they were testing method 1) or update
expected values.

Add NEWS.md entry: “Breaking-ish: default comb.method changed from 1 to
2 to match the paper’s recommendation.”

Bump DESCRIPTION patch (changes user-visible default).

CHECKPOINT 2B.i: review.

------------------------------------------------------------------------

### \[ \] 2C. Rename mislabeled `method_berger_boos`

1.  Investigate

Read `R/comparison_methods.R:539-563`. The function takes `pmax` of two
one-sided CIs at `alpha/2`. That is a Bonferroni union bound, not the
Berger-Boos profile-of-the-nuisance procedure (Berger and Boos 1994
maximize a tail probability over a confidence set for a nuisance
parameter; that is not what this function does).

Find every caller: Grep for `method_berger_boos` across
`/Users/jwbowers/repos/CMRSS_jake` and
`/Users/jwbowers/repos/combined_stephenson_tests`.

Check `tests/testthat/test-comparison_methods.R:77-115` (existing tests
that exercise this function).

Decide name with Jake: `method_bonferroni_two_sided` (clearest),
`method_min_p_two_sided`, or `method_pmax_two_sided`.

Decide whether to also implement an actual Berger-Boos nuisance
procedure as a new benchmark (likely separate task).

2.  Tests-first

In `tests/testthat/test-comparison_methods.R`, add:

- `test_that("renamed method exists and works")` that calls the new
  name.
- `test_that("old name method_berger_boos still works but warns deprecation")`
  that captures the warning message and verifies the result equals the
  new name.
- `test_that("renamed method computes element-wise pmax of two alpha/2 one-sided CIs")`
  to nail the actual semantics with a small synthetic example.

Run, observe failures. Pause.

CHECKPOINT 2C.t: Jake confirms the new name and decides on whether to
schedule a real Berger-Boos benchmark.

3.  Implementation

Rename the function in `R/comparison_methods.R:539`. Update its roxygen
`@title`, `@description`, `@references` (drop Berger and Boos 1994 from
the references for the renamed function), and `@examples`.

Add a thin wrapper
`method_berger_boos <- function(...) { .Deprecated("method_bonferroni_two_sided"); method_bonferroni_two_sided(...) }`
in the same file. Export the wrapper for one release cycle.

Update NAMESPACE (via `devtools::document()`).

Update test file header comment
(`tests/testthat/test-comparison_methods.R:1`).

Update any benchmark scripts in `combined_stephenson_tests/code/` (they
live in a separate repo; coordinate via a manual step or skip if Jake
handles separately).

Bump DESCRIPTION patch (renamed export).

Add NEWS.md entry.

CHECKPOINT 2C.i: review.

------------------------------------------------------------------------

## PRIORITY 3 – Cosmetic but worth doing

### \[ \] 3A. Comment the `j = m_s - k_s` relabeling in `solve_optimization`

1.  Investigate

Read `solve_optimization` and `HiGHS_sol_com` (`R/solvers.R:104-252`)
carefully. The variable index `j` in the LP corresponds to “number of
treated units in stratum s with effect \> c.” The paper indexes by
`k_s`. Identify the relabeling line(s) and confirm.

2.  Tests-first: none (pure comment); but verify no existing test relies
    on internal variable naming.

&nbsp;

3.  Implementation

Add a WHY block-comment near the constraint-construction loop in
`HiGHS_sol_com` and `Gurobi_sol_com` clarifying the index-name mapping
to the paper.

CHECKPOINT 3A: review (cheap; can batch with 3B and 3E).

### \[ \] 3B. Comment that LP coverage constraint is `<= p`, not `= p`

1.  Investigate

`HiGHS_sol_com` (`R/solvers.R:215`) uses `lhs = -Inf`, `rhs = p` for the
coverage row, i.e. `<=`. The paper’s eq:constraint_ks is written as
equality. The relaxation to `<=` is justified by monotonicity of the
objective in the number of “Inf-augmented” units.

2.  Tests-first: optional. Add a sanity test that the LP solution
    exhausts the budget when the optimum is attained at the boundary.

&nbsp;

3.  Implementation

WHY comment near `R/solvers.R:210-215` and the analogous Gurobi line
`R/solvers.R:364-365` (“=” vs “\<=”; the paper has equality, code uses
`<=` because the objective is monotone, citing the relevant proof step).

CHECKPOINT 3B: review.

### \[ \] 3C. Docstring fix for `set = "all"` alpha semantics

1.  Investigate

Read `com_conf_quant_larger_cre` (`R/CMRSS_CRE.R:725-797`, especially
line 773 and 791 where it splits `alpha/2`) and
`com_block_conf_quant_larger` (`R/CMRSS_SRE.R:1416-1485`, lines 1466 and
1479). Both pass `alpha = alpha/2` to the per-side helper, giving a
joint level-`alpha` CI (Bonferroni). The paper’s framing is
`1 - 2*alpha`. Document that the function’s `alpha` is the joint level
(1-alpha simultaneous coverage).

2.  Tests-first

Add `tests/testthat/test-set-all-coverage.R`:

- Synthetic null where the truth is known (e.g., sharp null,
  random-permutation truth).
- `test_that("set = 'all' alpha argument is the joint level")` verifies
  empirical coverage of the joint CI is at least `1 - alpha`.

3.  Implementation

Update the `@param alpha` blocks in both functions to make the
joint-level convention explicit. Note the conservativeness vs the
paper’s `1 - 2*alpha` framing.

CHECKPOINT 3C: review.

### \[ \] 3D. Harmonize defaults across CRE and SRE

1.  Investigate

`null.max`: SRE `pval_comb_block` defaults to `10^5`
(`R/CMRSS_SRE.R:1020`); SRE `com_block_conf_quant_larger` defaults to
`10^4` (`R/CMRSS_SRE.R:1424`); CRE `com_conf_quant_larger_cre` uses
`nperm = 10^4` (`R/CMRSS_CRE.R:726`).

`set`: CRE `com_conf_quant_larger_cre` defaults to `"treat"`
(`R/CMRSS_CRE.R:727`); SRE `com_block_conf_quant_larger` defaults to
`"all"` (`R/CMRSS_SRE.R:1418`).

Decide unified defaults with Jake. Recommend `null.max = 10^4` for
inversion functions (CI loops are expensive) and document; `set = "all"`
for both wrappers (more conservative, more useful default).

2.  Tests-first

Update tests that depend on the old defaults to pass them explicitly.

Add a `test_that("documented defaults are stable")` smoke test that
locks in the defaults via
[`formals()`](https://rdrr.io/r/base/formals.html).

3.  Implementation

Edit defaults in `R/CMRSS_SRE.R:1418, 1420, 1424` and
`R/CMRSS_CRE.R:726-727` to the chosen values.

Update roxygen `@param` text. NEWS.md entry. Bump DESCRIPTION patch
(default change).

CHECKPOINT 3D: review.

### \[ \] 3E. Docstring note for `Polynomial(std = FALSE)` vs paper eq:polynomial

1.  Investigate

Read `rank_score` in `R/CMRSS_CRE.R:84-110`. With `std = TRUE`, the
score is `(c(1:n)/(n+1))^(r-1)`. With `std = FALSE`, the score is
`(1:n)^(r-1)`, no Puri normalization.

Check `main.tex` eq:polynomial to see which convention the paper
assumes.

2.  Tests-first: not strictly necessary; a small test that
    `rank_score(n, list(name = "Polynomial", r = 2, std = TRUE))`
    matches the paper formula is helpful and cheap.

&nbsp;

3.  Implementation

Add a WHY comment at the relevant lines of `R/CMRSS_CRE.R`, plus an
`@details` paragraph in the roxygen for `rank_score` linking the paper
equation to each branch.

CHECKPOINT 3E: review.

### \[ \] 3F. Add missing test invariants identified by the audit

1.  Investigate

Confirm 1A and 1B tests already pin the `p` and `k` semantics via brute
force.

Confirm 2B test already asserts `comb.method = 2` is at least as
powerful as 1 in the Theorem 5 case.

Confirm 3C test already does the `set = "all"` coverage check.

2.  Tests-first / Implementation (this task is just bookkeeping if
    1A/1B/2B/3C tests are done)

If any of those three invariants are not yet covered, add the missing
tests in a single file `tests/testthat/test-paper-invariants.R`.
Otherwise, mark this item closed and reference the tests written in 1A,
1B, 2B, 3C.

CHECKPOINT 3F: review.

------------------------------------------------------------------------

## Cross-cutting build discipline (reapply at every implementation step)

Before each implementation phase, write tests first. Run them. Confirm
they fail (or pin current behavior) before writing code.

After each implementation, run `devtools::document()` then
`devtools::check()`.

Bump DESCRIPTION patch (currently 0.2.5) for: 1A (if numbers change), 1B
(if numbers change), 2A option (i) or (iii), 2B, 2C, 3D. Pure-comment
items (3A, 3B, 3E, 3F) do not require a bump.

Update `NEWS.md` for any user-visible default or signature change.

Pause for Jake review at every CHECKPOINT marker above.

------------------------------------------------------------------------

## Critical Files for Implementation

- `/Users/jwbowers/repos/CMRSS_jake/R/CMRSS_SRE.R`
- `/Users/jwbowers/repos/CMRSS_jake/R/CMRSS_CRE.R`
- `/Users/jwbowers/repos/CMRSS_jake/R/solvers.R`
- `/Users/jwbowers/repos/CMRSS_jake/R/comparison_methods.R`
- `/Users/jwbowers/repos/CMRSS_jake/tests/testthat/test-CMRSS_SRE.R`
