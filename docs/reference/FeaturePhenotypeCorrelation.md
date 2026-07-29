# FeaturePhenotypeCorrelation

This function performs correlation analysis of a specific feature
against a phenotype in the metadata. It can use linear mixed models
(LMM), simple linear regression (LM), or Pearson correlation. The
default regression line of the plot uses linear model

## Usage

``` r
FeaturePhenotypeCorrelation(
  mmo,
  feature_id,
  phenotype,
  groups,
  model = "lm",
  normalization = "Z",
  outdir = "FeaturePhenotypeCorrelation",
  width = 6,
  height = 6,
  save_output = TRUE
)
```

## Arguments

  - mmo:
    
    The `mmo` with feature data and metadata

  - feature\_id:
    
    The id of the feature to analyze

  - phenotype:
    
    The name of the phenotype in the metadata

  - groups:
    
    A vector of group names from the metadata containing phenotype data

  - model:
    
    The type of regression model to use. Options are 'lmm' for linear
    mixed model, 'lm' for simple linear regression, or 'pearson' for
    Pearson correlation (default: 'lmm')

  - normalization:
    
    The normalization method to use for feature data. Options are
    'None', 'Log', 'Meancentered', or 'Z' (default: 'Z')

  - outdir:
    
    The directory to save the output files

  - width:
    
    The width of the output plot in inches (default: 6)

  - height:
    
    The height of the output plot in inches (default: 6)

  - save\_output:
    
    A logical value indicating whether to save the output plot (default:
    TRUE)

## Value

A list of the plot and the raw data
