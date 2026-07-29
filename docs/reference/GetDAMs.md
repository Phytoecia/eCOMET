# Generates lists of DAMs (Differentially Accumulated Metabolites) for each comparison in the `mmo`

This function generates lists of upregulated and downregulated DAMs for
each pairwise comparison in the `mmo`. It uses log2 fold change and
adjusted p-value thresholds to determine significance. Make sure to run
PairwiseComp() for all desired comparisons before using this function.

## Usage

``` r
GetDAMs(mmo, fc_cutoff = 0.5849625, pval_cutoff = 0.05, use_padj = TRUE)
```

## Arguments

  - mmo:
    
    The `mmo` with pairwise comparison matrix

  - fc\_cutoff:
    
    The threshold of log2 fold change to be considered significant
    (default: 0.5849625, which is log2(1.5))

  - pval\_cutoff:
    
    The threshold of adjusted p-value to be considered significant
    (default: 0.05)

  - use\_padj:
    
    Boolean value indicating whether to use adjusted p-value (default:
    TRUE)

## Value

A list containing two lists: DAMs\_up and DAMs\_down

## Examples

``` r
if (FALSE) {
dams <- GetDAMs(mmo, fc_cutoff = 0.5849625, pval_cutoff = 0.05)
dams_up <- dams$DAMs_up
dams_down <- dams$DAMs_down
}
```
