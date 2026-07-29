# Use sample mass in the metadata file to normalize the peak area

This function normalizes the peak area in the feature data of the `mmo`
by the mass of each sample, provided in the metadata. The feature data
is replaced by (original value \* mean mass) / sample mass.

## Usage

``` r
MassNormalization(mmo)
```

## Arguments

  - mmo:
    
    The `mmo`

## Value

The `mmo` with normalized feature data (mmo$feature\_data)

## Examples

``` r
if (FALSE) {
mmo <- MassNormalization(mmo)
}
```
