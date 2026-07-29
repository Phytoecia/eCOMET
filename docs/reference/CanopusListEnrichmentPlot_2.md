# Generate a plot for enrichment analysis of Canopus-predicted terms across multiple levels

This function generates a plot for enrichment analysis of
Canopus-predicted terms across multiple levels, showing fold enrichment,
p-value, and subset count for each term level.

## Usage

``` r
CanopusListEnrichmentPlot_2(
  mmo,
  id_list,
  pthr = 0.05,
  outdir,
  height = 5,
  width = 5,
  topn = 5,
  pval = "pval",
  save_output = TRUE
)
```

## Arguments

  - mmo:
    
    The `mmo` with sirius annotation and normalized data

  - id\_list:
    
    A vector containing names of features to analyze

  - pthr:
    
    The threshold for adjusted p-value to be considered significant
    (default: 0.05)

  - outdir:
    
    The output file path for the enrichment plot

  - height:
    
    The height of the output plot in inches (default: 5)

  - width:
    
    The width of the output plot in inches (default: 5)

  - topn:
    
    The number of top terms to display in the plot (default: 5)

  - pval:
    
    pvalue options-pval or fdr (default: 'pval')

  - save\_output:
    
    boolean, whether to save the output plot (default: TRUE)

## Value

A list containing the enrichment plot and the enrichment results

## Examples

``` r
if (FALSE) {
CanopusListEnrichmentPlot_2(
 mmo, id_list = DAMs_up$control_vs_treatment1.up,
 pthr = 0.05, outdir = 'canopus_enrichment_plot_topn.pdf',
 height = 5, width = 5, topn = 5
)
}
```
