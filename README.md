# Combining Multiple Rank Sum Statistics

> R package for conducting randomization inference for quantiles of individual treatment effects, using combined rank sum statistics, both for completely randomized and stratified randomized experiments

## Installation

```
devtools::install_github("davidk91919/CMRSS")
```

## Load the Package

```
library(CMRSS)
```

## Explanation for Main Functions

```
?comb_p_val_cre ## get p-value from combined statistics in CRE
?com_conf_quant_larger_cre ## get the confidence interval for quantiles of individual effects, using combined statistics in CRE
?pval_comb_block ## get p-value from combined statistics in SRE
?com_block_conf_quant_larger ## get the confidence interval for quantiles of individual effects, using combined statistics in SRE
```
