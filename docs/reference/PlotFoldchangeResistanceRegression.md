# PlotFoldchangeResistanceRegression

This function plots the regression results of a feature against a fold
change in resistance, including regression line, p-value, and R-squared
value.

## Usage

``` r
PlotFoldchangeResistanceRegression(
  performance_regression,
  fold_change,
  color,
  outdir,
  width = 6,
  height = 6,
  save_output = TRUE
)
```

## Arguments

  - performance\_regression:
    
    The regression results data frame containing effect size, fold
    change, and tag. The output from GetPerformanceFeatureRegression,
    GetPerformanceFeatureLMM, or GetPerformanceFeatureCorrelation.

  - fold\_change:
    
    The name of the fold change column in the performance\_regression
    dataframe

  - color:
    
    A vector of colors for the points in the plot

  - outdir:
    
    The output file path for the regression plot

  - width:
    
    The width of the output plot in inches (default: 6)

  - height:
    
    The height of the output plot in inches (default: 6)

  - save\_output:
    
    A logical value indicating whether to save the output plot (default:
    TRUE)

## Value

A list containing the regression plot and the performance regression
data frame
