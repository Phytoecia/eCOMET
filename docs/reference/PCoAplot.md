# PCoAplot

This function generates a Principal Coordinates Analysis (PCoA) plot
based on the provided beta diversity distance matrix. It also performs
PERMANOVA analysis to assess the significance of group differences and
saves the results to CSV files.

## Usage

``` r
PCoAplot(
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
    
    The beta diversity distance matrix, output of GetBetaDiversity

  - outdir:
    
    The prefix for the output files

  - width:
    
    The width of the output PCoA plot (default: 6

  - height:
    
    The height of the output PCoA plot (default: 6)

  - color:
    
    A vector of colors for the points in the plot

  - save\_output:
    
    A logical value indicating whether to save the output plot (default:
    TRUE)

## Value

A list containing the PCoA plot, the PCoA coordinates, and the PERMANOVA
results

## Examples

``` r
if (FALSE) {
beta_diversity <- GetBetaDiversity(mmo, method = 'Gen.Uni',
 normalization = 'None', distance = 'dreams', filter_id = FALSE)
PCoAplot(mmo, betadiv = beta_diversity, outdir = 'output/PCoA', width = 6, height = 6)
}
```
