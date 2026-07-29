# Perform pairwise comparison between two groups in the `mmo`

This function performs pairwise comparison between two groups in the
`mmo`, calculating log2 fold change and adjusted p-values for given
comparison of two groups. The function adds the results to the
`mmo$pairwise` data frame.

## Usage

``` r
PairwiseComp(mmo, group1, group2, correction = "BH")
```

## Arguments

  - mmo:
    
    The `mmo`

  - group1:
    
    The name of the nominator group

  - group2:
    
    The name of the denominator group

  - correction:
    
    The method for multiple comparison correction. Options are 'BH',
    'holm', 'bonferroni', etc. Inherits from p.adjust() (default: 'BH')

## Value

The `mmo` with pairwise comparison results added to `mmo$pairwise`

## Examples

``` r
if (FALSE) {
mmo <- PairwiseComp(mmo, group1 = 'Control', group2 = 'Treatment1')
}
```
