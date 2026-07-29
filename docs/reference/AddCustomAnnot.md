# Add custom annotations to a `mmo`

Match features to a custom DB by m/z (ppm) and RT (minutes) tolerances
and attach a list-column of candidate compound IDs per feature.

## Usage

``` r
AddCustomAnnot(mmo, DB_file, mztol = 5, rttol = 0.5)
```

## Arguments

  - mmo:
    
    A `mmo` created by `GetMZmineFeature()`.

  - DB\_file:
    
    CSV path with at least columns `compound`, `mz`, `rt`.

  - mztol:
    
    m/z tolerance in ppm (default 5).

  - rttol:
    
    RT tolerance in minutes (default 0.5).

## Value

The same `mmo` with `mmo$custom_annot` (id, feature, custom\_annot).

## Examples

``` r
if (FALSE) {
mmo <- AddCustomAnnot(mmo, DB_file = "path/to/custom_db.csv", mztol = 5, rttol = 0.5)
}
```
