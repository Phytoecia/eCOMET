# Generate a plot for enrichment analysis of Canopus-predicted terms at a specific level using a list of vectors of features

This function generates a plot for enrichment analysis of
Canopus-predicted terms at a specific level, showing fold enrichment,
p-value, and subset count for each term.

## Usage

``` r
CanopusLevelEnrichmentPlot(
  mmo = mmo,
  comp.list,
  term_level = "NPC_pathway",
  pthr = 0.1,
  representation = "greater",
  outdir = "enrichment",
  height = 5,
  width = 5,
  pval = "pval",
  save_output = TRUE
)
```

## Arguments

  - mmo:
    
    The `mmo` with sirius annotation and normalized data

  - comp.list:
    
    A list to analyze, where each element is a vector of feature ids

  - term\_level:
    
    The level of term to use for enrichment analysis. Options are
    'NPC\_pathway', 'NPC\_superclass', 'NPC\_class',
    'ClassyFire\_superclass', 'ClassyFire\_class',
    'ClassyFire\_subclass', 'ClassyFire\_level5', or
    'ClassyFire\_most\_specific' (default: 'NPC\_pathway')

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
    
    The height of the output plot in inches (default: 5)

  - width:
    
    The width of the output plot in inches (default: 5)

  - pval:
    
    pvalue options-pval or fdr (default: 'pval')

  - save\_output:
    
    boolean, whether to save the output plot (default: TRUE)

## Value

A list containing the enrichment plot and the enrichment results

## Examples

``` r
if (FALSE) {
# Perform enrichment analysis for multiple comparisons using NPC_pathway level
comp.list <- list(
  comparison1 = DAMs_up$control_vs_treatment1.up,
  comparison2 = DAMs_up$control_vs_treatment2.up
)
CanopusLevelEnrichmentPlot(
 mmo, comp.list = comp.list, term_level = 'NPC_pathway',
 pthr = 0.1, representation = 'greater', outdir = 'enrichment_plot',
 height = 5, width = 5
)
}
```
