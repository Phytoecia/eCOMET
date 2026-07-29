# Convert feature IDs to names in the `mmo`

This function converts feature IDs to their corresponding names in the
`mmo`.

## Usage

``` r
IDToFeature(mmo, feature_ids)
```

## Arguments

  - mmo:
    
    The `mmo`

  - feature\_ids:
    
    A vector of feature IDs to convert

## Value

A vector of feature names corresponding to the input feature IDs

## Examples

``` r
if (FALSE) {
feature_names <- IDToFeature(mmo, feature_ids = c("1219", "2250", "3360"))
feature_names <- IDToFeature(mmo, feature_ids = mmo$feature_data$id[1:10])
feature_names <- IDToFeature(mmo,
 feature_ids = FeatureToID(mmo, feature_names = Glucosinolates)
) # if Glucosinolates is a vector of feature names
}
```
