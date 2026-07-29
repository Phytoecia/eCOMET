# Metabolite Set Enrichment Analysis (MSEA)

This function performs Metabolite Set Enrichment Analysis (MSEA) using
the fgsea package. It takes a ranked list of feature scores and tests
for enrichment of metabolite sets based on Canopus-predicted terms. The
results are saved as a CSV file and a PDF plot.

## Usage

``` r
MSEA(
  mmo,
  feature_id,
  feature_score,
  term_level = "NPC_class",
  pthr = 0.05,
  outdir = "MSEA",
  width = 8,
  height = 12,
  sig = FALSE,
  save_output = TRUE
)
```

## Arguments

  - mmo:
    
    The `mmo` with sirius annotation and normalized data

  - feature\_id:
    
    A vector of feature ids corresponding to the feature scores

  - feature\_score:
    
    A vector of feature scores (e.g., log2 fold changes)

  - term\_level:
    
    The level of term to use for enrichment analysis. Options are
    'NPC\_pathway', 'NPC\_superclass', 'NPC\_class',
    'ClassyFire\_superclass', 'ClassyFire\_class',
    'ClassyFire\_subclass', 'ClassyFire\_level5', or
    'ClassyFire\_most\_specific' (default: 'NPC\_class')

  - pthr:
    
    The threshold for adjusted p-value to be considered significant
    (default: 0.05)

  - outdir:
    
    The directory to save the output files (default: 'MSEA')

  - width:
    
    The width of the output plot in inches (default: 8)

  - height:
    
    The height of the output plot in inches (default: 12)

  - sig:
    
    A logical value indicating whether to return only significant terms
    (default: FALSE)

  - save\_output:
    
    A logical value indicating whether to save the output plot (default:
    TRUE)

## Value

A list of the plot and the enrichment results

## Examples

``` r
if (FALSE) {
# Perform MSEA using NPC_class level
MSEA(
 mmo, feature_name = rownames(DE_results), feature_score = DE_results$log2FoldChange,
 term_level = 'NPC_class', pthr = 0.05, outdir = 'MSEA_NPC_class',
 width = 8, height = 12, sig = FALSE
)
}
```
