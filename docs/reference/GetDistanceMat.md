# Get the distance matrix from the `mmo` based on the specified distance metric

Retrieve a feature distance matrix from the `mmo`

## Usage

``` r
GetDistanceMat(mmo, distance = "dreams")
```

## Arguments

  - mmo:
    
    The `mmo`

  - distance:
    
    Name of the distance matrix to retrieve. Built-in options are
    `'dreams'`, `'cosine'`, and `'m2ds'`. Any name passed to
    `AddCustomDist(mmo, name = ...)` is also valid.

## Value

The distance matrix (numeric matrix with feature ID row/col names).

## Details

Looks up a dissimilarity matrix stored in the `mmo` by name. Works with
the three built-in matrices added by `AddChemDist()` as well as any
custom matrix added via `AddCustomDist()`.

## Examples

``` r
if (FALSE) {
distance_matrix <- GetDistanceMat(mmo, distance = 'dreams')
distance_matrix <- GetDistanceMat(mmo, distance = 'tanimoto')  # custom
}
```
