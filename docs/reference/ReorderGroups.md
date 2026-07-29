# Reorder samples in the `mmo` based on group order

This function reorders the samples in the `mmo` based on a specified
group order. The function updates the order of samples in the feature
data, log-normalized data, z-score data, and mean-centered data. Use
this function before plotting heatmaps or other visualizations to ensure
consistent group ordering.

## Usage

``` r
ReorderGroups(mmo, group_order)
```

## Arguments

  - mmo:
    
    The `mmo`

  - group\_order:
    
    A vector specifying the desired order of groups

## Value

The `mmo` with reordered samples

## Examples

``` r
if (FALSE) {
mmo <- ReorderGroups(mmo, group_order = c("Control", "Treatment1", "Treatment2"))
}
```
