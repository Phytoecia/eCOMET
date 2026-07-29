# GetPerformanceFeatureCorrelation

This function calculates the Pearson correlation between each feature
and a specified phenotype performance in the metadata.

## Usage

``` r
GetPerformanceFeatureCorrelation(
  mmo,
  phenotype,
  groups,
  DAM.list,
  comparisons,
  cor_method = "pearson",
  normalization = "None"
)
```

## Arguments

  - mmo:
    
    The `mmo` with feature data and metadata

  - phenotype:
    
    The name of the phenotype performance in the metadata

  - groups:
    
    A vector of group names from the metadata containing performance
    data

  - DAM.list:
    
    A list of DAMs to tag features

  - comparisons:
    
    A list of pairwise comparisons to add fold change columns

  - cor\_method:
    
    The correlation method to use. Options are 'pearson', 'spearman', or
    'kendall' (default: 'pearson')

  - normalization:
    
    The normalization method to use for feature data. Options are
    'None', 'Log', 'Meancentered', or 'Z' (default: 'None')

## Value

A data frame containing correlation results for each feature, including
effect size (correlation coefficient), p-value, and fold change columns
for specified comparisons.
