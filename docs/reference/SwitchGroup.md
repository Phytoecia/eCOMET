# Switch the group column in the `mmo`

This function switches the group column in the metadata of the `mmo` to
a new specified column. The new group column must exist in the metadata
file.

## Usage

``` r
SwitchGroup(mmo, new_group_col)
```

## Arguments

  - mmo:
    
    The `mmo`

  - new\_group\_col:
    
    The name of the new group column in the metadata file

## Value

The `mmo` with the updated group column

## Examples

``` r
if (FALSE) {
mmo <- SwitchGroup(mmo, new_group_col = "genotype")
}
```
