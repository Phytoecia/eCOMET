# Enrichment analysis for Canopus-predicted terms

This function performs enrichment analysis for Canopus-predicted terms
on a given list of features.

## Usage

``` r
CanopusLevelEnrichmentAnal(
  mmo,
  list_test,
  pthr = 0.1,
  sig = TRUE,
  term_level = "NPC_pathway",
  representation = "greater",
  pval = "pval"
)
```

## Arguments

  - mmo:
    
    The `mmo` with sirius annotation and normalized data

  - list\_test:
    
    A vector containing ids of features to analyze

  - pthr:
    
    The threshold for adjusted p-value to be considered significant
    (default: 0.1)

  - sig:
    
    A logical value indicating whether to return only significant terms
    (default: TRUE)

  - term\_level:
    
    The level of term to use for enrichment analysis Options are
    'NPC\_pathway', 'NPC\_superclass', 'NPC\_class',
    'ClassyFire\_superclass', 'ClassyFire\_class',
    'ClassyFire\_subclass', 'ClassyFire\_level5', or
    'ClassyFire\_most\_specific' (default: 'NPC\_pathway')

  - representation:
    
    The representation type for enrichment analysis. Options are
    'greater' for overrepresentation (default: 'greater')

  - pval:
    
    pvalue options-pval or fdr (default: 'pval')

## Value

A data frame containing the enrichment results, including term level,
term name, subset count, total count, fold enrichment, p-value, and
adjusted p-value (FDR)

## Examples

``` r
if (FALSE) {
# Perform enrichment analysis for a list of features using NPC_pathway level
sig_terms <- CanopusLevelEnrichmentAnal(
 mmo, list_test = c("feature1", "feature2"), pthr = 0.1,
 sig = TRUE, term_level = 'NPC_pathway', representation = 'greater'
)
# Perform enrichment analysis for a list of features using
# ClassyFire_class level and return all terms
all_terms <- CanopusLevelEnrichmentAnal(
 mmo, list_test = c("feature1", "feature2"), pthr = 0.1,
 sig = FALSE, term_level = 'ClassyFire_class', representation = 'greater'
)
}
```
