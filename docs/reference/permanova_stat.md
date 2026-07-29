# Perform PERMANOVA and pairwise comparisons

This function performs PERMANOVA on the given data and metadata, with
options for filtering groups. It also conducts post-hoc pairwise
comparisons and adjusts p-values for multiple testing. The function
returns the PERMANOVA results, raw pairwise comparison results, and
matrices of adjusted p-values, F values, and R square for pairwise
comparisons

## Usage

``` r
permanova_stat(
  data,
  metadata,
  mode,
  filter_group = FALSE,
  group_list = NULL,
  permutations = 5000
)
```

## Arguments

  - data:
    
    A data frame or distance matrix for PERMANOVA

  - metadata:
    
    A data frame containing sample metadata, including a 'group' column

  - mode:
    
    The mode of the input data: 'data' for raw data or 'distance' for a
    distance matrix

  - filter\_group:
    
    Boolean to filter groups based on a provided list (default: FALSE)

  - group\_list:
    
    A vector of group names to filter (default: NULL)

  - permutations:
    
    The number of permutations for the PERMANOVA test (default: 5000)

## Value

A list containing the PERMANOVA results, raw pairwise comparison
results, and matrices of adjusted p-values, F values, and R square for
pairwise comparisons

## Examples

``` r
if (FALSE) {
permanova_results <- permanova_stat(
 data = feature_data, metadata = mmo$metadata,
 mode = 'data', filter_group = TRUE, group_list = c("Control", "Treatment1"),
 permutations = 5000
)
permanova_results <- permanova_stat(
 data = betadiv, metadata = mmo$metadata,
 mode = 'distance', permutations = 10000
)
}
```
