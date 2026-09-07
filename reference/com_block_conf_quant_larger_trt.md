# Helper function for Simultaneous Inference for multiple quantiles on SRE

Output is a lower limit of prediction intervals for prespecified
quantiles under the all-units convention. Used internally by
[`com_block_conf_quant_larger`](https://bowers-illinois-edu.github.io/CMRSS/reference/com_block_conf_quant_larger.md).

## Usage

``` r
com_block_conf_quant_larger_trt(
  Z,
  Y,
  block,
  k.vec = NULL,
  methods.list.all,
  weight.name = "asymp.opt",
  opt.method = "ILP_auto",
  comb.method = 1,
  stat.null = NULL,
  null.max = 10^4,
  Z.perm = NULL,
  tol = 0.01,
  alpha = 0.1
)
```

## Arguments

- Z:

  An \\n\\ dimensional treatment assignment vector.

- Y:

  An \\n\\ dimensional observed outcome vector.

- block:

  An \\n\\ dimensional vector specifying block of each units.

- k.vec:

  Optional vector of specific quantile indices to compute. If NULL,
  computes for all quantiles.

- methods.list.all:

  A list of lists of lists. Corresponds to the method for each stratum
  for each different stratified rank sum statistic.

- weight.name:

  Weighting method to be implemented. If "asymp.opt", asymptotically
  optimal scheme under a class of local alternatives is adjusted, if
  "dist.free", design-free scheme is adjusted.

- opt.method:

  Optimization method. Options include:

  - "ILP_gurobi": Integer linear programming with Gurobi solver

  - "LP_gurobi": Linear programming relaxation with Gurobi solver

  - "ILP_highs": Integer linear programming with HiGHS solver
    (open-source)

  - "LP_highs": Linear programming relaxation with HiGHS solver

  - "ILP" or "ILP_auto": Integer LP with auto-selected solver

  - "LP" or "LP_auto": Linear programming with auto-selected solver

  HiGHS is recommended for users without a Gurobi license. Both solvers
  produce equivalent results.

- comb.method:

  Integer specifying the combination method:

  - 1: Combine statistics across strata first, then take maximum across
    methods (default)

  - 2: Take maximum across methods within each stratum first, then
    combine across strata

- stat.null:

  An vector whose empirical distribution approximates the randomization
  distribution of the combined stratified rank sum statistic.

- null.max:

  A positive integer representing the number of permutations for
  approximating the randomization distribution of the rank sum
  statistic.

- tol:

  Tolerance for root-finding algorithm.

- alpha:

  Significance level for confidence intervals.

## Value

Vector of length `n` giving lower confidence limits for the all-units
quantiles \\\tau\_{(1)}, \ldots, \tau\_{(n)}\\. Only entries indexed by
\\k = n - m, \ldots, n\\ (where \\m = \mathrm{sum}(Z)\\) are computed by
the LP; lower entries are padded.

## Details

The hypothesis tested at each `k` is the all-units null \\H\_{k, c}:
\tau\_{(k)} \le c\\, where \\\tau\_{(1)} \le \tau\_{(2)} \le \ldots \le
\tau\_{(n)}\\ are the sorted individual treatment effects across all
\\n\\ units. The LP coverage constraint sets `p <- n - k` (lines 1185
and 1248 of `R/CMRSS_SRE.R`); `k` ranges over \\1, \ldots, n\\ and the
loop fills `quantiles[k]` for \\k = n - m, \ldots, n\\, where \\m =
\mathrm{sum}(Z)\\, padding the remaining lower entries with
`quantiles[n - m]`.

By the validity result extending ZL24quantile (cf. main.tex eq:H_kc_t
and the discussion at lines 479 and 519-520), the inverted \\1 -
\alpha\\ lower bound for \\\tau\_{(k)}\\ is simultaneously a valid \\1 -
\alpha\\ lower prediction bound for \\\tau\_{\mathrm{treat}(k -
n\_{\mathrm{control}})}\\, where \\n\_{\mathrm{control}} = n - m\\.
Equivalently, the LP minimum at all-units index \\k_R\\ matches the LP
minimum that
[`pval_comb_block`](https://bowers-illinois-edu.github.io/CMRSS/reference/pval_comb_block.md)
(treated-only convention, `p <- m - k`) would compute at \\k_C = k_R -
n\_{\mathrm{control}}\\. The two functions therefore produce the same
numerical confidence bounds despite indexing their `k` differently.

[`com_block_conf_quant_larger`](https://bowers-illinois-edu.github.io/CMRSS/reference/com_block_conf_quant_larger.md)
calls this function in three modes: directly with `(Z, Y)` for
`set = "treat"`, and with `(1 - Z, -Y)` for `set = "control"` and the
control half of `set = "all"`. Under the swap, individual effects are
preserved via the relabeling \\\tilde{\tau}\_{\mathrm{treat}(k)} =
\tau\_{\mathrm{control}(k)}\\ (cf. main.tex line 479), so the swap-space
lower CIs returned by this function are directly lower CIs on
\\\tau\_{\mathrm{control}(k)}\\ in the original problem – no sign flip
or index reversal is needed.
