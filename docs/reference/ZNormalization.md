# Z-normalize the peak area in the `mmo`

This function applies Z-score normalization to the peak area in the
feature data of the `mmo`. Z-scores are calculated per feature (row)
across samples. Features with zero variance cannot be Z-normalized and
are returned as NA.

## Usage

``` r
ZNormalization(mmo, imputed_data = FALSE)
```

## Arguments

  - mmo:
    
    The `mmo`

  - imputed\_data:
    
    Whether to use imputed feature data (default = FALSE)

## Value

The `mmo` with Z-normalized feature data stored in `mmo$zscore`
