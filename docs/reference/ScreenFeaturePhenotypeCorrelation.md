# Screen feature-phenotype correlation

Use metadata-provided variables (any phenotypes or environmental
variables) to screen feature-phenotype correlation linear model, linear
mixed model (using groups as random effect), or correlation (Pearson,
Spearman, Kendall) are supported

## Usage

``` r
ScreenFeaturePhenotypeCorrelation(
  mmo,
  phenotype,
  groups,
  model = "lm",
  normalization = "None"
)
```

## Arguments

  - mmo:
    
    The `mmo` with feature data and metadata

  - phenotype:
    
    The name of the phenotype in the metadata

  - groups:
    
    A vector of group names from the metadata containing phenotype data

  - model:
    
    The type of regression model to use. Options are 'lmm' for linear
    mixed model, 'lm' for simple linear regression, or 'pearson',
    'spearman', 'kendall' for correlation (default: 'lm')

  - normalization:
    
    The normalization method to use for feature data. Options are
    'None', 'Log', 'Meancentered', or 'Z' (default: 'Z')

## Value

A list of the plot and the raw data
