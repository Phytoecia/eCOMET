# Convert feature abundances to presence / absence

This function converts the feature abundance matrix in a `mmo` into a
binary presence/absence matrix and stores it as a new component of the
`mmo` (`mmo$feature_presence`).

## Usage

``` r
FeaturePresence(mmo, threshold = 1)
```

## Arguments

  - mmo:
    
    The `mmo`

  - threshold:
    
    Numeric threshold for presence (default = 1). Values \> threshold
    are set to 1, values \<= threshold or NA are set to 0.

## Value

The `mmo` with a new presence/absence table stored in
`mmo$feature_presence`

## Details

A feature is considered present (1) if its abundance is greater than a
specified threshold, and absent (0) otherwise.

This function does NOT overwrite `mmo$feature_data`.

## Examples

``` r
if (FALSE) {
mmo <- FeaturePresence(mmo, threshold = 1)
}
```
