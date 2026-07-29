# Calculate group means from the `mmo`

This function calculates and returns a dataframe of mean feature values
for each group in the `mmo`, with options for normalization and
filtering. Use SwitchGroup() to change the grouping variable before
running this function.

## Usage

``` r
GetGroupMeans(
  mmo,
  normalization = "None",
  filter_id = FALSE,
  id_list = NULL,
  filter_group = FALSE,
  group_list = NULL
)
```

## Arguments

  - mmo:
    
    The `mmo`

  - normalization:
    
    The normalization method to use. Options are 'None', 'Log',
    'Meancentered', or 'Z'

  - filter\_id:
    
    Boolean to filter features based on a provided list (default: FALSE)

  - id\_list:
    
    A vector of feature names to filter (default: NULL)

  - filter\_group:
    
    Boolean to filter groups based on a provided list (default: FALSE)

  - group\_list:
    
    A vector of group names to filter (default: NULL)

## Value

A data frame containing the mean feature values for each group

## Examples

``` r
if (FALSE) {
group_means <- GetGroupMeans(mmo, normalization = 'Log')
group_means <- GetGroupMeans(mmo,
 normalization = 'None',
 filter_id = TRUE, id_list = Glucosinolates
) # if Glucosinolates is a vector of feature names
group_means <- GetGroupMeans(mmo,
 normalization = 'Z',
 filter_group = TRUE,
 group_list = c("Control", "Treatment1")
)
group_means <- GetGroupMeans(mmo,
 normalization = 'Meancentered',
 filter_id = TRUE, id_list = Glucosinolates,
 filter_group = TRUE, group_list = c("Control", "Treatment1")
) # if Glucosinolates is a vector of feature names
}
```
