# Calculate log2 fold change for a given control group

This function calculates and returns a dataframe of log2 fold change
values for each group compared to a specified control group. Takes
inputs from GetGroupMeans() function.

## Usage

``` r
GetLog2FoldChange(group_means, control_group)
```

## Arguments

  - group\_means:
    
    A data frame containing the mean feature values for each group

  - control\_group:
    
    The name of the control group to compare against

## Value

A data frame with log2 fold change values for each group compared to the
control group

## Examples

``` r
if (FALSE) {
fold_change <- GetLog2FoldChange(GetGroupMeans(mmo,
                 normalization = 'Log'),
                 control_group = 'Control')
}
```
