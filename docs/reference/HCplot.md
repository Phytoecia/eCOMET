# HCplot

Hierarchical clustering of samples from a precomputed beta-diversity
matrix, plotted as a phylogram with tip labels colored by
`mmo$metadata$group`.

## Usage

``` r
HCplot(
  mmo,
  betadiv,
  outdir,
  hclust_method = "average",
  color = NULL,
  tip_label_size = 2.5,
  width = 10,
  height = 7,
  save_output = TRUE
)
```

## Arguments

  - mmo:
    
    The `mmo`

  - betadiv:
    
    Beta diversity distance matrix, output of GetBetaDiversity()

  - outdir:
    
    Output file prefix (e.g., "output/HC"). PNG saved here if
    save\_output=TRUE.

  - hclust\_method:
    
    hclust linkage method (default: "average"; alternatives:
    "complete","ward.D2")

  - color:
    
    Named vector of colors for groups. If NULL, uses Set3 palette.

  - tip\_label\_size:
    
    Tip label text size (default: 2.5)

  - width:
    
    Plot width in inches (default: 10)

  - height:
    
    Plot height in inches (default: 7)

  - save\_output:
    
    Whether to save the plot (default: TRUE)

## Value

A list containing: plot (ggplot object), hc (hclust), phy (phylo)

## Examples

``` r
if (FALSE) {
bet <- GetBetaDiversity(mmo, method = 'bray', normalization = 'None')
hc <- HCplot(mmo, betadiv = bet, outdir = "output/HC_bray")
hc$plot
}
```
