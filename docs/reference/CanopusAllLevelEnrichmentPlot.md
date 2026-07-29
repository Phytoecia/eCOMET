# Generate a plot for enrichment analysis of Canopus-predicted terms across all levels

This function generates a plot for enrichment analysis of
Canopus-predicted terms across all levels, showing fold enrichment,
p-value, and subset count for each term level.

## Usage

``` r
CanopusAllLevelEnrichmentPlot(
  mmo = mmo,
  comp.list,
  terms = "all_terms",
  term_levels = NULL,
  pthr = 0.1,
  representation = "greater",
  outdir = "enrichment",
  height = 10,
  width = 8,
  pval = "pval",
  save_output = TRUE
)
```

## Arguments

  - mmo:
    
    The `mmo` with sirius annotation and normalized data

  - comp.list:
    
    A list to analyze, where each element is a vector of feature names

  - terms:
    
    The terms to analyze. Options are 'all\_terms', 'NPC', 'ClassyFire',
    or 'custom' (default: 'all\_terms')

  - term\_levels:
    
    list of custom term levels to use

  - pthr:
    
    The threshold for adjusted p-value to be considered significant
    (default: 0.1)

  - representation:
    
    The representation type for enrichment analysis. Options are
    'greater' for overrepresentation (default: 'greater')

  - outdir:
    
    The output directory for saving the plot and the enrichment results
    (default: 'enrichment')

  - height:
    
    The height of the output plot in inches (default: 10)

  - width:
    
    The width of the output plot in inches (default: 8)

  - pval:
    
    pvalue options-pval or fdr (default: 'pval')

  - save\_output:
    
    boolean, whether to save the output plot (default: TRUE)

## Value

A list of the plot and the enrichment results

## Examples

``` r
if (FALSE) {
comp.list <- list(
  comparison1 = DAMs_up$control_vs_treatment1.up,
  comparison2 = DAMs_up$control_vs_treatment2.up
)
CanopusAllLevelEnrichmentPlot(
 mmo, comp.list = comp.list, terms = 'all_terms',
 pthr = 0.1, representation = 'greater', outdir = 'enrichment_all_levels',
 height = 10, width = 8
)
CanopusAllLevelEnrichmentPlot(
 mmo, comp.list = comp.list, terms = 'NPC',
 pthr = 0.1, representation = 'greater', outdir = 'enrichment_NPC_levels',
 height = 10, width = 8
)
CanopusAllLevelEnrichmentPlot(
 mmo, comp.list = comp.list, terms = 'ClassyFire',
 pthr = 0.1, representation = 'greater', outdir = 'enrichment_ClassyFire_levels',
 height = 10, width = 8
)
}
```
