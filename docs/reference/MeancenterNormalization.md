# Mean-center the peak area in the `mmo`

This function applies mean-centering to the peak area in the feature
data of the `mmo`. Mean-centering is performed per feature (row) across
samples. Features with zero variance are returned as all zeros and are
reported in a warning.

## Usage

``` r
MeancenterNormalization(mmo, imputed_data = FALSE)
```

## Arguments

  - mmo:
    
    The `mmo`

  - imputed\_data:
    
    Whether to use imputed feature data (default = FALSE)

## Value

The `mmo` with mean-centered feature data stored in `mmo$meancentered`
