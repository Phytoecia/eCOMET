# Volcano plot for visualizing differential metabolite analysis results

This function generates a volcano plot using data from `mmo$pairwise`
(`PairwiseComp(mmo, 'group1', 'group2')` should be specified),
highlighting upregulated and downregulated features based on log2 fold
change and adjusted p-value

## Usage

``` r
VolcanoPlot(
  mmo,
  comp,
  topk = 10,
  log2FC_thr = 1,
  pthr = 0.05,
  outdir = "volcano.png",
  height = 5,
  width = 5,
  save_output = TRUE,
  use_padj = TRUE
)
```

## Arguments

  - mmo:
    
    The `mmo` with pairwise comparison matrix

  - comp:
    
    The comparison to visualize, e.g., 'group1\_vs\_group2

  - topk:
    
    The number of top features to label in the plot (default: 10)

  - log2FC\_thr:
    
    The threshold of log2 fold change to be considered significant
    (default: 1)

  - pthr:
    
    The threshold of adjusted p-value to be considered significant
    (default: 0.05)

  - outdir:
    
    The output file path for the volcano plot (default: 'volcano.png')

  - height:
    
    The height of the output plot in inches (default: 5)

  - width:
    
    The width of the output plot in inches (default: 5)

  - save\_output:
    
    A logical value indicating whether to save the output plot (default:
    TRUE)

  - use\_padj:
    
    A logical value indicating whether to use adjusted p-value (default:
    TRUE)

## Value

A list containing the volcano plot and the data used to generate it

## Examples

``` r
if (FALSE) {
VolcanoPlot(
 mmo, comp = 'Control_vs_Treatment1',
 topk = 10, log2FC_thr = 1, pthr = 0.05,
 outdir = 'volcano_con_tre1.png', height = 5, width = 5
)
}
```
