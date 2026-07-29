# \#' Replace zero and NA values in the `mmo`

This function replaces zero values in the feature table of a `mmo`.
Imputed data are stored in mmo$imputed\_feature\_data, to be used for
downsteam analyses including PairwiseComp(). Note that imputation
affects normalizations (Log-transformation, etc.), as well as chemical
diversity calculations that uses presence/absence.

## Usage

``` r
ReplaceZero(mmo, method = c("one", "half_min"))
```

## Arguments

  - mmo:
    
    A `mmo` containing `feature_data`

  - method:
    
    Replacement method:
    
      - "one": replace zeros and NA values with 1
    
      - "half\_min": replace zeros and NA values with half of the
        smallest non-zero value in each feature (row)

## Value

The updated `mmo`
