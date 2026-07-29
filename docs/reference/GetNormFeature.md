# Retrieve feature data from the `mmo`, with normalization options

This function retrieves the feature data from the `mmo` based on the
specified normalization method.

## Usage

``` r
GetNormFeature(mmo, normalization)
```

## Arguments

  - mmo:
    
    The `mmo`

  - normalization:
    
    The normalization method to use. Options are 'None','PA', 'Log',
    'Meancentered', 'Z', or 'Imputed'

## Value

The feature data corresponding to the specified normalization method

## Examples

``` r
if (FALSE) {
feature_data <- GetNormFeature(mmo, normalization = 'Log')
}
```
