# NMDSplot

This function generates a Non-metric Multidimensional Scaling (NMDS)
plot based on the provided beta diversity distance matrix. It also
performs PERMANOVA analysis to assess the significance of group
differences and saves the results to CSV files.

## Usage

``` r
NMDSplot(
  mmo,
  betadiv,
  outdir,
  width = 6,
  height = 6,
  color = NULL,
  save_output = TRUE
)
```

## Arguments

  - mmo:
    
    The `mmo` containing metadata

  - betadiv:
    
    The beta diversity distance matrix, output of GetBetaDiversity()

  - outdir:
    
    The outdir for the output files

  - width:
    
    The width of the output NMDS plot (default: 6)

  - height:
    
    The height of the output NMDS plot (default: 6)

  - color:
    
    A vector of colors for the groups in the plot

  - save\_output:
    
    A logical value indicating whether to save the output plot (default:
    TRUE)

## Value

A list containing the NMDS plot, the NMDS coordinates, and the PERMANOVA
results

## Examples

``` r
if (FALSE) {
beta_diversity <- GetBetaDiversity(mmo, method = 'Gen.Uni',
 normalization = 'None', distance = 'dreams', filter_id = FALSE)
# Use method = 'bray' or 'jaccard' if you want to use just feature abundance
# without considering feature spectral dissimilarity
NMDSplot(mmo, betadiv = beta_diversity, outdir = 'output/NMDS', width = 6, height = 6)
}
```
