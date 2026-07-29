# CalculateGroupBetaDistance

This function calculates the beta diversity distance between a reference
group and other groups in the metadata.

## Usage

``` r
CalculateGroupBetaDistance(mmo, beta_div, reference_group, groups)
```

## Arguments

  - mmo:
    
    The `mmo` containing feature data and metadata

  - beta\_div:
    
    The beta diversity distance matrix, output of GetBetaDiversity()

  - reference\_group:
    
    The name of the reference group to compare against

  - groups:
    
    A vector of group names from the metadata to calculate distances for

## Value

A data frame containing the group names, sample names, and their
corresponding beta diversity distances from the reference group.

## Examples

``` r
if (FALSE) {
beta_diversity <- GetBetaDiversity(mmo, method = 'Gen.Uni',
 normalization = 'None', distance = 'dreams', filter_id = FALSE)
group_distances <- CalculateGroupBetaDistance(mmo, beta_div = beta_diversity,
 reference_group = 'Control', groups = c('Control', 'Treatment1', 'Treatment2'))
}
```
