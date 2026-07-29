# GetPerformanceFeatureRegression

This function performs linear regression analysis of all features
against a phenotype performance in the metadata.

## Usage

``` r
GetPerformanceFeatureRegression(
  mmo,
  phenotype,
  groups,
  DAM.list,
  comparisons,
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

  - normalization:
    
    The normalization method to use for feature data. Options are
    'None', 'Log', 'Meancentered', or 'Z' (default: 'None')

## Value

A data frame containing regression results for each feature, including
effect size, p-value, and fold change columns for specified comparisons.
