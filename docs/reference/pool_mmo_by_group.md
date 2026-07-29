# pool\_mmo\_by\_group

Pools sample columns within each group into one pseudo-sample per group.
Feature rows are preserved (no filtering of features). Keeps all other
slots of the `mmo` unchanged by copying the `mmo` first.

## Usage

``` r
pool_mmo_by_group(mmo, group_col = "group")
```

## Arguments

  - mmo:
    
    A `mmo`

  - group\_col:
    
    column in `mmo$metadata` used for grouping (default: "group")

## Value

`mmo` with feature\_data containing one column per group and metadata
updated accordingly
