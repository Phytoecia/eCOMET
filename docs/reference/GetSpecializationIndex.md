# GetSpecializationIndex

This function calculates the specialization index for a given `mmo`,
normalization method, and optional filtering by groups and features.

## Usage

``` r
GetSpecializationIndex(
  mmo,
  normalization = "None",
  filter_group = FALSE,
  group_list = NULL,
  filter_id = FALSE,
  id_list = NULL
)
```

## Arguments

  - mmo:
    
    The `mmo` containing feature data and metadata

  - normalization:
    
    The normalization method to use for feature data. Options are
    'None', 'Log', 'Meancentered', or 'Z' (default: 'None')

  - filter\_group:
    
    A boolean indicating whether to filter the feature data by a
    specific group list (default: FALSE)

  - group\_list:
    
    A list of groups to filter the feature data by, if filter\_group is
    TRUE (default: NULL)

  - filter\_id:
    
    A boolean indicating whether to filter the feature data by a
    specific list (default: FALSE)

  - id\_list:
    
    A list of feature names to filter the feature data by, if filter\_id
    is TRUE (default: NULL)

## Value

A data frame containing the specialization index for each group in the
metadata, with columns for group and specialization index.

## Examples

``` r
if (FALSE) {
specialization_index <- GetSpecializationIndex(mmo,
                                               normalization = 'None',
                                               filter_group = FALSE)
specialization_index <- GetSpecializationIndex(mmo,
                                               normalization = 'Z',
                                               filter_group = TRUE,
                                               group_list = c('Control', 'Treatment1'),
                                               filter_id = TRUE)
}
```
