# Convert feature names to IDs in the `mmo`

This function converts feature names to their corresponding IDs in the
`mmo`.

## Usage

``` r
FeatureToID(mmo, feature_names)
```

## Arguments

  - mmo:
    
    The `mmo`

  - feature\_names:
    
    A vector of feature names to convert

## Value

A vector of feature IDs corresponding to the input feature names

## Examples

``` r
if (FALSE) {
feature_ids <- FeatureToID(mmo, feature_names = c("100.0_5.0", "150.0_10.0"))
feature_ids <- FeatureToID(mmo, feature_names = mmo$feature_data$feature[1:10])
feature_ids <- FeatureToID(mmo,
 feature_names = Glucosinolates
) # if Glucosinolates is a vector of feature names
}
```
