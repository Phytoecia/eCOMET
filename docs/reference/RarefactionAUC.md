# RarefactionAUC

This function calculates the rarefaction AUC for a given rarefied
richness table

## Usage

``` r
RarefactionAUC(rarefied_richness, n_boot = 1000, seed = 513)
```

## Arguments

  - rarefied\_richness:
    
    The rarefied richness object containing feature data, annotations,
    and pairwise comparisons

  - n\_boot:
    
    The number of bootstrap samples to use (default: 1000)

  - seed:
    
    The seed to use for the bootstrap samples (default: 513)

## Value

A list containing the rarefaction AUC for each group in the metadata,
with columns for group and AUC

## Examples

``` r
if (FALSE) {
RarefactionAUC(rarefied_richness, n_boot = 1000)
}
```
