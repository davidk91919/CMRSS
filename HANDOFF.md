# HANDOFF: CMRSS paper-vs-code audit

Originally written 2026-04-26 after a paper-vs-code audit. Updated
2026-04-28 to reflect a repo move, David Kim’s reply on item 1A, and a
follow-on finding that complicates David’s reply. Updated 2026-05-12 to
record CI infrastructure fixes and explicit skip markers for the 10
tests that needed to be revisited when David replied. Updated 2026-05-19
to record the closure of item 1A (treated-only convention adopted; tests
re-enabled in CMRSS 0.2.7) and scenario-(a) step 3
(`com_block_conf_quant_larger` wrapper audit; no bug; function
documentation added). Updated 2026-09-07 to record the merge with
`upstream/main`, the arrival of the revise-and-resubmit decision on the
paper, David Kim’s step back from the project, and a verification that
the paper still reproduces from its pinned CMRSS commit. Read top to
bottom.

## TL;DR for the next session (2026-09-07 onward)

Read this section first. The 2026-05-20 TL;DR below it is still accurate
about item 1A, but its `main` SHA and its paper-status paragraph have
both been overtaken; both are corrected in place.

### Where the repositories stand

`main` is `b198562` and is identical on the local clone, on `origin`
(`bowers-illinois-edu/CMRSS`), and on `upstream` (`davidk91919/CMRSS`).
The divergence that existed at the start of the session is gone: 13
ahead, 0 behind became 0 and 0.

`b198562` is a merge commit. Its second parent is David’s `b67c935`
(2026-05-25), which changed one line of `tests/testthat/test-pval-cre.R`
from `k <- floor(0.8 * n)` to `k <- floor(0.8 * m)`. The conflict was
resolved in favor of our version of that file, and the merge changed no
files: `git diff aa67001 b198562` is empty. The merge commit message
carries the derivation of why our version is the correct one; read it
before reopening the question.

Why David’s edit was not taken: it fixes the CMRSS side of that test but
leaves the RIQITE side alone, and RIQITE implements the all-units null
`H_{k,c}` of `main.tex` eq:H_kc (line 510, `k in 1..n`) while CMRSS
implements the treated-only null `H_{k,c}^treat` of eq:H_kc_t (line 459,
`k in 1..n_treat`). Passing one `k` to both compares different
hypotheses. At `n = 50`, `m = 25`, `k = 20`, RIQITE’s constraint on the
treated units is vacuous (`20 - 25 = -5`) and it returns exactly
1.00000, against CMRSS’s 0.77923 (Wilcoxon) and 0.71548 (Stephenson);
the test asserts agreement at `tolerance = 0.01`, so David’s version
fails by 0.22 and 0.28.

Jake has `push: true` on `davidk91919/CMRSS` and pushed `main` there
directly on 2026-09-07, so the canonical repo went from 0.2.5 to 0.2.7
in one step. He drafted a Slack note to David explaining the one place
they differ. He declined a dual push URL on `origin`, preferring to push
to the canonical repo deliberately, so `main` can drift again. Check
`git rev-list --count upstream/main..main` before assuming the two
agree.

`gh-pages` differs between the two repos and always will, since each
repo’s pkgdown workflow builds its own site. Leave it alone.

### The paper is at revise-and-resubmit

This supersedes the “manuscript under review” risk posture stated in the
2026-05-20 TL;DR and in the next-session menu at the bottom of this
file. The wait for the R&R signal is over.

The constraint that replaces it, in Jake’s words: he needs to reproduce
the results as submitted, exactly, before making any revision. So the
gating question for an item is no longer “does this change numerical
outputs before the R&R,” it is “has Jake reproduced and recorded the
submitted numbers yet.” Ask him before starting any item under “Possible
numerical change” in the menu below.

One operational rule follows, and it belongs to the paper repo rather
than this one: do not run
[`renv::snapshot()`](https://rstudio.github.io/renv/reference/snapshot.html)
or
[`renv::update()`](https://rstudio.github.io/renv/reference/update.html)
in `~/repos/combined_stephenson_tests` until the submitted numbers are
reproduced and saved. Either command would move the CMRSS pin off
`1371f315` and destroy the reference point.

### Replication is safe, and was verified rather than assumed

`~/repos/combined_stephenson_tests/renv.lock` pins CMRSS 0.2.5 at commit
`1371f315` from `bowers-illinois-edu/CMRSS` (2025-12-24). Confirmed
2026-09-07 via the GitHub API that the commit is still fetchable, so the
April archive-and-recreate of that repo did not break the pin.

More than that, the pin barely matters. Between `1371f315` and `b198562`
the only executable change anywhere in `R/` is the six-line range check
on `k` inside `pval_comb_block`; everything else is roxygen text plus
the `DESCRIPTION` bump. The paper’s four run scripts call
`com_block_conf_quant_larger` and `com_conf_quant_larger_cre`, which
reach the solver through `com_block_conf_quant_larger_trt` and never
route through `pval_comb_block`. So the guard cannot fire on the paper’s
code path, and current `main` reproduces the submitted numbers as well
as the pin does. This closes scenario-(a) step 5 (the paper-script
audit) as verified rather than deferred.

`code/sre_simulation_run.R` calls
[`library(gurobi)`](https://rdrr.io/r/base/library.html), so reproducing
the SRE simulation needs a live Gurobi license. The electric-teachers
scripts and the CRE simulation do not load it.

### David Kim is stepping back

He has taken an industry job, and the active authors are now Jake and
Xinran Li. Jake’s decision on 2026-09-07: `davidk91919/CMRSS` stays the
canonical home for now, so the repo-layout section further down this
file is still correct. Do not propose archiving or renaming without
asking.

### State of the tree

`devtools::test()` on `b198562` reports 0 failures, 0 warnings, 150
tests passing, and 4 skipped, all four being the Gurobi-gated
comparisons in `test-solvers.R` that skip because Gurobi is not
installed on this machine. `devtools::check()` was not re-run this
session because the merge changed no files and `aa67001` was already
check-clean at 0 errors, 0 warnings, and 2 repo-hygiene NOTEs.

### What was not done

Nothing in `R/` or `tests/` changed this session. The open PLAN items
are exactly as the menu at the bottom of this file lists them, with the
R&R gating rewritten. Scenario-(a) step 4 (refactoring
`com_block_conf_quant_larger_trt` from all-units `k` to treated-only
`k`) is still deferred and still Jake’s call.

## TL;DR for the next session (2026-05-20 onward)

**Item 1A is closed.** Treated-only convention adopted; 10 previously
skipped tests are unskipped and passing; the cross-validation against
RIQITE uses the `k_R = k_C + (n - m)` translation under a shared
`Z.perm`. `com_block_conf_quant_larger`’s `set = "all"` and
`set = "control"` paths were audited against `main.tex` and found
correct; tests and function documentation are in place.
`devtools::check()` is clean (0 errors, 0 warnings, 2 pre-existing
repo-hygiene NOTEs). `main` was at `99af2eb` when this paragraph was
written; as of 2026-09-07 it is `b198562`. Working tree should be clean.

**Status of the older “paused indefinitely” framing below: stale.** The
2026-04-28 section and several below it describe a state of waiting for
David Kim’s reply on item 1A. That wait is over and the resolution is
recorded in the 2026-05-19 update. Read those older sections only for
forensic context, not for next-step guidance.

Paper status, as written on 2026-05-19 and now superseded: the
manuscript was under review, and the risk posture was to avoid
numerical-output changes until the R&R signal arrived. It has arrived.
See “The paper is at revise-and-resubmit” in the 2026-09-07 section
above for the constraint that replaces it.

**The “Suggested next session opening” section at the bottom of this
file is up to date** (rewritten 2026-05-19). Start there.

## Update 2026-05-19: item 1A closed (treated-only convention; 0.2.7)

The “paused indefinitely” status described below is over for item 1A.
Jake confirmed the resolution: `R/CMRSS_SRE.R:1034` (`p <- m - k`) is
correct, `pval_comb_block` tests the treated-only hypothesis
`H_{k,c}^treat` with `k in 1..sum(Z)`, and CMRSS and RIQITE must return
the same p-values when called consistently. The 10 skipped tests have
been re-enabled.

What was changed in 0.2.7 (uncommitted at the time of this update; see
`git status` for the staged set):

1.  `R/CMRSS_SRE.R`: `pval_comb_block` roxygen and example rewritten to
    document the treated-only convention. Stale `# p <- N - k` comment
    removed.
2.  `man/pval_comb_block.Rd`: regenerated by `devtools::document()`.
3.  `DESCRIPTION`: version 0.2.6 -\> 0.2.7.
4.  `NEWS.md`: 0.2.7 entry describing the documentation and test
    changes.
5.  10 test sites previously marked
    `skip("Pending k-convention resolution; ...")` had the skip removed
    and their `k` formulas converted from the all-units convention
    (`floor(0.8 * N)`, `floor(0.9 * n)`, etc.) to the treated-only
    convention (`floor(0.8 * s * m)`, `floor(0.9 * sum(Z))`, etc.).
6.  4 additional gurobi-gated sites in `test-solvers.R` (lines 119, 257,
    367, 439) had the same latent k-convention bug; their `k` formulas
    were updated in the same pass so they will pass on a machine with
    Gurobi.
7.  The two RIQITE cross-validation tests (`test-pval-cre.R` and the
    second `test_that` in `test-pval_scre.R`) were rewritten to:
    - Pre-compute a shared `Z.perm` matrix and pass it to both packages
      (via
      [`RIQITE::null_dist`](https://rdrr.io/pkg/RIQITE/man/null_dist.html)
      and `com_null_dist_block`).
    - Apply the LP-equivalence translation
      `k_R (RIQITE all-units) = k_C (CMRSS treated-only) + (n - m)`.
    - Assert agreement at `tolerance = 1e-6`.

    Reason the translation is needed: `H_{k_R,c}` (RIQITE all-units) and
    `H_{k_C,c}^treat` (CMRSS treated-only) are different hypotheses, but
    their LP minima for the rank-sum statistic coincide under that
    translation because controls’ xi does not enter `Z*xi` in the test
    statistic. Freely assigning all `(n - m)` controls to `xi = c`
    absorbs `(n - m)` of RIQITE’s all-units constraint slack, leaving
    exactly `k_C = k_R - (n - m)` treated units constrained. The test
    verifies this with `diff = 0` across multiple `k` values when the
    permutation set is shared.
8.  `tests/testthat/test-pval-comb-block-p-convention.R` header trimmed:
    the “inversion sibling uses `p <- n - k`” note was removed; the
    `com_block_conf_quant_larger_trt` reconciliation is now tracked only
    in PLAN.md and the “Still open” list below.

`devtools::check()` on the result: 0 errors, 0 warnings, 2 NOTEs. Both
NOTEs are pre-existing repo-hygiene items (hidden files like `.lintr`,
`.claude`; non-standard top-level files like `HANDOFF.md`, `AGENTS.md`,
`CLAUDE.md`, `PLAN.md`, `Makefile`, `load.R`). `devtools::test()`: 135+
pass, 4 skip (all gurobi-gated; no k-convention skips remain).

### Scenario-(a) step 3 closed 2026-05-19 (no functional bug)

The `com_block_conf_quant_larger` `set = "all"` and `set = "control"`
paths were audited against `main.tex` (eq:H_kc_t, eq:H_kc_c, paper lines
459-466, 479, 502, 940, 948-950). Finding: the wrapper is mathematically
correct.

- Under the swap `Z <- 1 - Z, Y <- -Y`, individual effects are preserved
  via the relabeling `tilde_tau_treat(k) = tau_control(k)` (paper line
  479). The swap-space lower CIs returned by
  `com_block_conf_quant_larger_trt(1 - Z, -Y, ...)` are *directly* lower
  CIs on `tau_control(k)` in the original problem – no sign flip or
  index reversal is needed. My earlier suspicion that the wrapper was
  returning “negated upper bounds” was wrong; I was thinking of the swap
  as `tilde_tau = -tau`, but the correct potential-outcome bookkeeping
  (`tilde_Y(1) = -Y(0)`, `tilde_Y(0) = -Y(1)`) leaves `tilde_tau = tau`
  unchanged.
- `set = "all"` is the Bonferroni pool from paper line 502: it calls the
  treat and control branches at `alpha/2` each, sorts, and returns. The
  user’s `alpha` is the joint coverage level (1 - alpha simultaneous
  coverage), consistent with paper line 950.

Work landed this pass:

- New file `tests/testthat/test-com_block_conf_quant_larger-sets.R` with
  4 `test_that` blocks (9 assertions). Pins length and type, pins
  `set = "control"` equals direct `_trt(1 - Z, -Y, ...)` with the
  correct slice, pins `set = "all"` equals
  `sort(c(treat_branch, control_branch))` at `alpha/2` each, pins the
  Bonferroni alpha split. All checks use `set.seed` paired calls so the
  RNG consumed for the null permutations is identical between the
  wrapper call and the manual replication; equality is exact
  (`tolerance = 1e-6`), not MC-tolerant.
- Function documentation on `com_block_conf_quant_larger_trt` expanded
  with `@details` explaining (i) the all-units convention
  (`p <- n - k`), (ii) the validity-theorem equivalence to
  `pval_comb_block`’s treated-only convention under the translation
  `k_R = k_C + (n - m)`, (iii) the swap semantics for the control / all
  branches. Internal function only, no exported-symbol change, no
  version bump.

### Still open from the larger scenario-(a) audit

1.  **`com_block_conf_quant_larger_trt` at `R/CMRSS_SRE.R:1171, 1234`
    still uses `p <- n - k` (all-units).** Now correctly *documented* as
    the all-units convention with the LP-equivalence note. But the
    asymmetry with `pval_comb_block` (treated-only) remains a
    maintenance hazard. Optional refactor: change `_trt` to take
    treated-only `k` internally and update the wrapper’s slicing
    accordingly. Functionally equivalent; risk is that any downstream
    caller of `_trt` (which is `@keywords internal` but still callable)
    would break.
2.  **Paper experiment scripts in
    `/Users/jwbowers/repos/combined_stephenson_tests/code/`.** Per the
    original audit, the scripts call only the wrapper
    `com_block_conf_quant_larger` and `com_conf_quant_larger_cre`, not
    `pval_comb_block` directly. With the wrapper now confirmed correct,
    paper-side numbers are unaffected. Audit is low priority; defer.

Other items still in PLAN.md (1B, 2A, 2B, 2C, 3A-3F) are unchanged by
this session. With item 1A closed, item 2B (default `comb.method` flip
from 1 to 2) is now unblocked.

## Update 2026-05-12: CI now runs on push (10 tests skipped pending David)

The substantive pause described below is unchanged. What changed today
is that GitHub Actions is now wired up properly on
`bowers-illinois-edu/CMRSS`, and the 10 tests blocked by item 1A carry
explicit skip markers so future-Jake can find them with a single grep
when the conversation with David resumes.

What was broken: push events on `main` were not triggering any
workflows. The pkgdown site had been deployed by hand from
`davidk91919/CMRSS` to `gh-pages` rather than by CI. Root cause was a
UI-only “I understand my workflows, go ahead and enable them” gate that
had never been clicked at the repo Actions page. The API reported
workflows as `state: active` regardless, so this was not visible from
`gh`.

What was changed:

1.  Jake clicked the enable gate in the browser; the org Actions policy
    was confirmed as “Allow all actions and reusable workflows.”
2.  `workflow_dispatch:` added to `.github/workflows/R-CMD-check.yaml`
    so manual triggers work (parity with pkgdown.yaml). Commit
    `a3b62dd`.
3.  `any::highs` added to `extra-packages:` in R-CMD-check.yaml. Without
    it, `dependencies: '"hard"'` skipped highs (in Suggests) and
    [`get_default_solver()`](https://bowers-illinois-edu.github.io/CMRSS/reference/get_default_solver.md)’s
    example aborted with “No solver available.” Commit `ebebfa5`.
4.  Ten tests `skip()`’d with the message
    `"Pending k-convention resolution; see PLAN.md item 1A and commit 53001b0"`.
    These are the same ten test sites flagged in commit `53001b0`’s
    message; their `pval_comb_block` calls use the “all-units k =
    floor(0.8 \* N)” convention, which the guard at `R/CMRSS_SRE.R:1034`
    rejects. Commit `108745d`.

Skipped tests, by file and `test_that` header line:

- `tests/testthat/test-CMRSS_SRE.R` (5): lines 136, 182, 319, 352,
  407. 
- `tests/testthat/test-pval_scre.R` (2): lines 6, 58.
- `tests/testthat/test-pval-cre.R` (1): line 8.
- `tests/testthat/test-solvers.R` (2): lines 57, 199.

Current matrix result: 5/5 OS configs pass (macos-latest,
windows-latest, ubuntu-latest devel/release/oldrel-1). The pkgdown
workflow rebuilds the site on every push to `main` and deploys to
`gh-pages`, so manual
[`pkgdown::deploy_to_branch()`](https://pkgdown.r-lib.org/reference/deploy_to_branch.html)
is no longer required.

When David replies and item 1A is resolved, find these tests with:

    grep -rn "Pending k-convention resolution" tests/testthat/

Remove the `skip()` line at the top of each affected `test_that` block.
If the resolution is “tests are wrong, the guard is right,” also adjust
each test’s `k = floor(0.8 * N)` to a value within `1..sum(Z)`. If the
resolution is “guard is wrong, tests are right,” remove the guard at
`R/CMRSS_SRE.R:1034` and the dedicated tests in
`tests/testthat/test-pval-comb-block-p-convention.R`.

## Current state (2026-04-28, still applies)

Work is **paused indefinitely**, awaiting David Kim’s reply on a
follow-up question. Jake’s framing in the Slack draft is “tackle this
after the dissertation is submitted,” so the pause may run weeks or
months. A fresh Claude picking this up should not assume the pause is
short.

The two commits from this session are pushed to
`bowers-illinois-edu/CMRSS` (a fork of `davidk91919/CMRSS`):

    53001b0 Validate k in pval_comb_block; record David's reply on item 1A
    a3d09fa Add renv-based dependency management for R 4.6.0

The Slack follow-up message to David is drafted at
`/tmp/slack_to_david_followup.md` (not in the repo). Jake has edited it;
the version on disk is the one to send. As of this write Jake had not
yet sent it. When the next session starts, check whether the file still
exists, whether David has replied (Slack is the channel, not email), and
update the file or this handoff accordingly.

## Repo layout (post 2026-04-28 move)

The package’s canonical home is now `davidk91919/CMRSS`. Prior to
2026-04-28 it lived at `bowers-illinois-edu/CMRSS`; that org repo was
mirror-pushed to David’s account, renamed to
`bowers-illinois-edu/CMRSS_archive`, and archived (read-only). An
`archive` tag and `archive` branch on the archive repo mark the pre-move
HEAD `f2703fb` for posterity.

`bowers-illinois-edu/CMRSS` was then re-created as a fresh GitHub fork
of `davidk91919/CMRSS`. The local clone at
`/Users/jwbowers/repos/CMRSS_jake` has:

- `origin` -\> `git@github.com:bowers-illinois-edu/CMRSS.git` (the org
  fork)
- `upstream` -\> `git@github.com:davidk91919/CMRSS.git` (David’s
  canonical)
- `main` tracks `origin/main`

Standard fork workflow: branch off `main`, push to `origin`, open PRs
against `upstream/main`. Pull David’s work with
`git fetch upstream && git merge upstream/main`.

## The open question for David (item 1A, second pass)

David replied earlier on 2026-04-28 that line `R/CMRSS_SRE.R:1034`
(`p <- m - k`, introduced in his commit `762c4d08` on 2025-12-18) is
correct as code. He said the bug is in the documentation – the docstring
claiming “k between 1 and n” – and that he would fix the docs upstream.
Item 1A in `PLAN.md` was marked resolved at that point.

Subsequent investigation in the same session turned up a stronger signal
that complicates David’s reply:

- `tests/testthat/test-pval-cre.R:43-56` cross-validates
  [`CMRSS::pval_comb_block`](https://bowers-illinois-edu.github.io/CMRSS/reference/pval_comb_block.md)
  against
  [`RIQITE::pval_quantile`](https://rdrr.io/pkg/RIQITE/man/pval_quantile.html)
  at the same `(k, c)` with `k = floor(0.8 * n)`. RIQITE uses the
  all-units convention. Pre-`762c4d08` this test passed. Post-, the call
  has `k = 40 > sum(Z) = 25`, the LP is infeasible, and the function
  silently returns `p.value = 0` – the cross-validation no longer
  cross-validates.
- The same all-units `k = floor(0.8 * N)` (or `floor(0.9 * n)`) pattern
  appears in **ten** existing test sites:
  - `tests/testthat/test-CMRSS_SRE.R:171, 216, 344, 386, 442`
  - `tests/testthat/test-pval_scre.R:42, 90`
  - `tests/testthat/test-pval-cre.R:49`
  - `tests/testthat/test-solvers.R:95, 233`

That is not a documentation problem. It is the function’s behavior
having shifted from `H_{k,c}` (all units) to `H_{k,c}^treat` (treated
only). The Slack follow-up to David asks him to choose between:

1.  **Intentional shift to treated-only.** The function now tests
    `H_{k,c}^treat`. Then we need to update the docstring, examples, the
    ten test sites, and audit anything downstream – in particular
    `com_block_conf_quant_larger`’s `set = "all"` path and the paper’s
    experiment scripts in
    `/Users/jwbowers/repos/combined_stephenson_tests/code/` – for
    assumptions about all-units `k`. (Per the original 2026-04-26 audit,
    the paper’s experiment scripts call only the wrapper
    `com_block_conf_quant_larger` and `com_conf_quant_larger_cre`, not
    `pval_comb_block` directly, so paper-side audit may be quick.)

2.  **`762c4d08` is a regression.** Line 1034 should revert to
    `p <- N - k`. The test suite goes back to green and there is no
    downstream audit.

**Both options are consistent with David’s “the docs are wrong” reply in
some reading**, which is why the follow-up is needed before more code
moves.

## What was done in this session that survives

All committed and pushed to `bowers-illinois-edu/CMRSS`:

1.  **renv setup** (commit `a3d09fa`). R 4.6.0 project library,
    `renv.lock`, `.Rprofile`, `renv/activate.R`, `renv/settings.json`,
    `renv/.gitignore`. `.Rbuildignore` updated to exclude renv from
    package build. ~305 packages installed; `gurobi` excluded
    (commercial, not on CRAN).
    [`renv::status()`](https://rstudio.github.io/renv/reference/status.html)
    reports clean.

2.  **Input validation in `pval_comb_block`** (commit `53001b0`). Added
    near `R/CMRSS_SRE.R:1036` (just below the `p <- m - k` line):

    ``` r

    if (!is.numeric(k) || length(k) != 1L || k < 1 || k > m) {
      stop(sprintf(
        "k = %s is outside the feasible range 1..sum(Z) = %d for the current LP coverage constraint p = m - k at R/CMRSS_SRE.R:1034.",
        format(k), m
      ))
    }
    ```

    The error message is deliberately neutral about the hypothesis – it
    states the constraint without claiming `H_{k,c}` or `H_{k,c}^treat`.
    Correct under either resolution of the David question. Bumps
    `DESCRIPTION` to 0.2.6 and adds a `NEWS.md` entry.

3.  **Test rewrite at
    `tests/testthat/test-pval-comb-block-p-convention.R`** (commit
    `53001b0`). Replaces the original three skipped all-units tests with
    three executable treated-only invariants:

    1.  no silent `p.value = 0 + test.stat = Inf` on out-of-range k
        (passes once the input validation guard is in place);
    2.  valid k in `1..m_total` yields finite stat and p in `[0, 1]`;
    3.  boundary `k = m_total` (so `p = 0`) is feasible. All three pass
        under the current code with the input-validation guard. They are
        neutral between (a) and (b).

4.  **HANDOFF.md and PLAN.md** updated to record David’s first reply and
    (in this rewrite) the follow-on finding.

Test edits I had **deliberately reverted** before committing:
`tests/testthat/test-CMRSS_SRE.R:344-345, 386-388, 442-444` (changing
`floor(0.8 * N)` to `floor(0.8 * m_total)`). These were assuming
interpretation (a). Reverting makes the failure pattern uniform across
all ten sites, which is easier for David to reproduce.

## How David should reproduce

The Slack message tells him:

``` bash
git clone git@github.com:bowers-illinois-edu/CMRSS.git CMRSS_jake_audit
cd CMRSS_jake_audit
Rscript -e "renv::restore()"
Rscript -e "devtools::load_all('.', quiet = TRUE); testthat::test_dir('tests/testthat')"
```

He will see ten errors, all of the same shape:

    Error in pval_comb_block(...):
    k = 40 is outside the feasible range 1..sum(Z) = 25 for the current LP
    coverage constraint p = m - k at R/CMRSS_SRE.R:1034.

The new test file `test-pval-comb-block-p-convention.R` passes 9/9 –
neutral signal, just confirms the guard works.

## Restarting after David replies

Likely scenario (a) – intentional treated-only shift. Steps:

1.  Update `pval_comb_block` docstring and example to use treated-only
    `k` (k in `1..sum(Z)`, example `k = floor(0.9 * sum(Z))` or
    similar). Match wording to David’s preferred framing
    (`H_{k,c}^treat`).
2.  Update the ten test sites listed above. Most just need
    `floor(0.8 * N)` -\> `floor(0.8 * sum(Z))` or `floor(0.9 * n)` -\>
    `floor(0.9 * sum(Z))`.
3.  **Audit `com_block_conf_quant_larger`’s `set = "all"` path.** The
    wrapper at `R/CMRSS_SRE.R:~1422` swaps `Z <- 1 - Z` for the
    “control” and “all” branches, then calls into the same machinery.
    After the swap, `m` and `n - m` switch roles. Whether the wrapper
    correctly translates the user’s intended (all-units) `k` into a
    treated-only `k` for the post-swap call needs to be checked
    carefully – this is the most likely place for a subtle bug to hide.
4.  **Decide what to do about `com_block_conf_quant_larger_trt`** at
    `R/CMRSS_SRE.R:1171, 1234`. Both lines still have `p <- n - k`
    (all-units). David did not touch them in `762c4d08`. If
    `pval_comb_block` is treated-only and the inversion sibling is
    all-units, they are testing different hypotheses despite the `_trt`
    suffix. Likely needs reconciliation.
5.  **Audit paper experiment scripts in
    `/Users/jwbowers/repos/combined_stephenson_tests/code/`** for any
    call to `com_block_conf_quant_larger` that might depend on the
    pre-2025-12-18 behavior. The paper’s reference code at
    `code/codes_20251026.R:1426` still has `p = N - k`. Re-running would
    probably be needed if the wrapper’s behavior changed.
6.  Update `tests/testthat/test-pval-comb-block-p-convention.R`’s header
    to remove “the inversion sibling uses `p <- n - k`” line if step 4
    changes it.
7.  `devtools::document()`, `devtools::check()`, and run the full test
    suite; expect 0 failures.

If scenario (b) – regression – steps:

1.  Revert `R/CMRSS_SRE.R:1034` to `p <- N - k`. Keep the input
    validation guard but update its error message and feasibility range
    to `k <= n` (matching `n - k >= 0`).
2.  The ten test sites and the new `test-pval-comb-block-p-convention.R`
    all pass without further edit. Update
    `test-pval-comb-block-p-convention.R` header to match.
3.  `devtools::document()`, `devtools::check()`.

If David replies with something Jake hasn’t anticipated: stop and check
with Jake, do not improvise.

## Other items still in PLAN.md (not touched this session)

Priority 1: - **1B**. Verify column range `0:nb[i]` in
`comb_matrix_block_stratum` (`R/CMRSS_SRE.R:704`) against paper’s
`eq:comb_per_stratum`. Possibly an over-enumeration. Independent of 1A;
can move on this if 1A is blocked.

Priority 2: - **2A**. Tie handling. `ties.method = "first"` in
[`rank()`](https://rdrr.io/r/base/rank.html) calls at
`R/CMRSS_CRE.R:13, 278` and `R/CMRSS_SRE.R:231, 799`. Paper may want
`"random"` for finite-sample validity with discrete outcomes
(`electric_teachers` has many ties). Fully independent of 1A. - **2B**.
Default `comb.method = 1` should flip to `2` per paper’s Theorem 5.
Depends on 1A being settled (the `p`/`k` semantics need to be fixed
first). - **2C**. `method_berger_boos` (`R/comparison_methods.R:539`) is
mislabeled; it does a Bonferroni `pmax`, not Berger-Boos. Rename.

Priority 3: - **3A**, **3B**, **3C**, **3D**, **3E**, **3F**. Cosmetic /
docstring fixes. See `PLAN.md` for detail.

## Important context to preserve

- **Jake’s role**: paper’s first author and an applied statistician at
  UIUC. Treat the paper at
  `/Users/jwbowers/repos/combined_stephenson_tests/main.tex` as the
  authoritative description of intended behavior. Hypotheses `H_{k,c}`
  (eq:H_kc, ~line 511) and `H_{k,c}^treat` (Section 1.2) are both
  legitimate; which one this function tests is the open question.
- **David Kim’s role**: co-author, current `cre` of the package, and
  active contributor. Do not assume his commits are wrong without
  evidence. Ask before reverting his work.
- **Paper’s experiment scripts** in
  `/Users/jwbowers/repos/combined_stephenson_tests/code/` (excluding the
  embedded `codes_20251026.R` and `old_codes/`) produced the published
  numbers. They call `com_block_conf_quant_larger` (the wrapper) and
  `com_conf_quant_larger_cre`, not `pval_comb_block` directly. Per the
  original audit, paper-side numbers are likely unaffected by
  `762c4d08`, but if scenario (a) wins, the wrapper’s `set = "all"`
  translation needs verification before claiming that.
- **Solvers**: HiGHS (open-source) is the default and what `R CMD check`
  will use unless Gurobi is installed locally. Tests use
  `skip_if_not(solver_available(...))`. Gurobi is in `Suggests` but is
  not in the renv lockfile.
- **renv**: project library is at
  `~/Library/Caches/org.R-project.R/R/renv/library/CMRSS_jake-2747cf55/...`
  (renv default cache layout). To restore on a fresh machine:
  `Rscript -e "renv::restore()"`. R 4.6.0 is the locked version.
- **ASCII only**: no unicode em dashes, en dashes, arrows, fancy quotes,
  ellipses, decorative bullets. Use `---`, `--`, `->`, straight quotes,
  `...`. Jake’s global rule.
- **Workflow rules** at
  `/Users/jwbowers/repos/ai_workflow/CLAUDE_CODING.md`: tests-first;
  pause for review at three checkpoints (after tests, after
  implementation, at any unresolved design decision); WHY comments only;
  bump `DESCRIPTION` patch when an exported symbol or user-visible
  default changes. Project memory at
  `/Users/jwbowers/.claude/projects/-Users-jwbowers-repos-CMRSS-jake/memory/feedback_coding_workflow.md`
  points to the same file.
- **Bowers-illinois-edu org**: Jake’s affiliation. The org repo is the
  working fork; David’s repo is canonical.

## Files outside the repo to remember

- `/tmp/slack_to_david_followup.md` – the active Slack draft to send.
  Edited by Jake; do not regenerate without checking with him.
- `/tmp/message_to_david.md`, `/tmp/slack_to_david.md` – earlier drafts
  from 2026-04-26 (pre-David-reply). Probably stale; check before using.
- `/Users/jwbowers/repos/combined_stephenson_tests/` – the paper repo.
  Main entry: `main.tex`. Reference R code: `code/codes_20251026.R`.
  Experiment scripts: rest of `code/`.

## Suggested next session opening (rewritten 2026-05-19; supersedes the older one)

1.  Read the **TL;DR at the top of this file**. Skim the 2026-05-19
    update. The older sections (2026-04-28 “Current state”, “The open
    question for David”, “Restarting after David replies”, etc.) are
    historical only – do not act on them.
2.  Confirm the working tree is clean and `main` is at the latest pushed
    commit (`git status`, `git log --oneline -3`).
3.  Re-read workflow memory at
    `/Users/jwbowers/.claude/projects/-Users-jwbowers-repos-CMRSS-jake/memory/feedback_coding_workflow.md`
    before any code work. Re-read the terminology memory at
    `feedback_terminology_function_documentation.md` (write “function
    documentation,” not “docstring”).
4.  Ask Jake which open item to tackle. The menu below was ordered by
    risk against an under-review manuscript. The paper is now at R&R, so
    the ordering still holds but the gate has changed: before starting
    anything under “Possible numerical change,” ask whether he has
    reproduced and recorded the submitted numbers yet.
    - **Safe (no numerical change):**
      - **2C** – rename `method_berger_boos` to e.g.
        `method_bonferroni_two_sided` (the function does Bonferroni
        `pmax`, not Berger-Boos; `R/comparison_methods.R:539-563`). Pure
        rename + deprecation shim.
      - **Coverage simulation for the wrapper** – verify the (1 - alpha)
        joint coverage claim on synthetic data with known truth. No code
        change, just a long test.
      - **3A, 3B, 3E** – pure function-documentation clarifications
        (PLAN.md).
      - **Scenario-(a) step 4** – optional refactor of
        `com_block_conf_quant_larger_trt` from all-units to treated-only
        `k`. Functionally equivalent; mechanically removes the asymmetry
        with `pval_comb_block`.
    - **Possible numerical change (gate: submitted numbers reproduced
      and recorded first):**
      - **1B** – column range `0:nb[i]` in `comb_matrix_block_stratum`
        (`R/CMRSS_SRE.R:704`) vs paper `eq:comb_per_stratum`. Possible
        over-enumeration. If wrong, `comb.method = 2` outputs shift.
      - **2A** – `ties.method = "first"` vs `"random"`
        (`R/CMRSS_CRE.R:13, 278`; `R/CMRSS_SRE.R:231, 799`).
        Finite-sample validity question; `electric_teachers` has many
        ties.
      - **2B** – default `comb.method = 1` -\> `2` (now unblocked by
        item 1A). User-visible default change.
      - **3C, 3D** – alpha-semantics docs and default harmonization that
        touch user-facing behavior.
    - **Scenario-(a) step 5** (paper-script audit at
      `combined_stephenson_tests/code/`) – done 2026-09-07. Verified
      that all four run scripts call `com_block_conf_quant_larger` or
      `com_conf_quant_larger_cre` and never `pval_comb_block`, so the
      `k` guard is unreachable from the paper’s code path. See the
      2026-09-07 section at the top of this file.
5.  Do not start coding before Jake picks an item and confirms scope.
    For any item that might change numerical outputs, ask first whether
    the submitted numbers have been reproduced and recorded.
