# MMO
mmo is desinged to get input from mzmine, sirius, and qemistry to perform downstream statistic analyses.


## 0. Input
Following input files would be imported
1) mzmine_feature.csv
CSV export of mzmine
2) canopus_formula_summary.tsv
3) structure_identification
4) metadata.csv
  metadata should have three or more columns; the first three being sample, group, and mass. Additional information can be provided
5) fingerprint_distance.csv
  output of qemistree 

## 1. Initializing the object
```
#Give directories
mzmine_featuredir <- "240913_mzmine_feature.csv"
metadatadir <- "metadata_test.csv"
canopus_formuladir <- "canopus_formula_summary.tsv"
canopus_structuredir <- "structure_identifications.tsv"
```
First import the mzmine feature list and metadata to create a `mmo` object.

`mmo <- GetMZmineFeature(mzmine_featuredir, metadatadir)`

then add annotation generated from sirius to the object

`mmo <- AddSiriusAnnot(mmo, canopus_structuredir, canopus_formuladir)`

then peform imputation

`mmo <- ReplaceZero(mmo) # Replace 0 and NA values by half of the minimum value`

then normalize using sample mass used for extraction

`mmo <- MassNormalization(mmo) # Normalize peak area by sample mass in metadata`

then generate various normalized through mean-centering, log, and z-score.
```
mmo <- MeancenterNormalization(mmo) # Add mean-centered area
mmo <- LogNormalization(mmo) # Add log-transformed area
mmo <- ZNormalization(mmo) # Add Zscore
```

### mmo object
The `mmo` object contains following information
1) Feature quantification value ($feature_data)
2) Metadata ($metadata)
3) Pairwise comparison data ($pairwise)
4) Log, mean-centered, and Z-transformed feature quantification value ($log, $meancentered, and #zscore)
5) Sirius annotation ($sirius_annot)

## 2. Summarise using PLS-DA
By using PLS-DA plot, the general distribution of metabolic profile of each sample and group can be visualized. 
For each group give colors

`custom_colors <- c("con" = "#999999", "sl1" = "#fdcdac", "px1" = "#cbd5e8", "px2" = "#86aaed", 
                   "px3" = "#3869c7", "sl2" = "#e89e6b", "sl3" = "#b65c1f", "le1" = "#b3e2cd", 
                   "le2" = "#67daa7", "le3" = "#1fb270")`

Then plot PLS-DA

`PLSDAplot(mmo, color = custom_colors, topk = 0,outdir = 'PLSDA_Z_withoutLoading.pdf', normalization = 'Z', filter = FALSE, filter_list = NULL)`

The `topk` parameter can be adjusted to plot loadings of top features on the plot. The `normalization` can be `None`, `Log`, `Meancentered`, and `Z`.
If a subset of features to be used for ploting, set `filter = TRUE` then provide the list of features to the `filter_list`

## 3. Pairwise Comparison (DAMs)
Many analyses targets to find **Differentially Accumulated Metabolites (DAMs; corresponding to the DEGs in transcriptomics)**. DAMs can be defined by thresholds of log2-fold change and adjusted p-value. Those two metrics are calculated by following code. Note that the divisor group is at the left.

`mmo <- PairwiseComp(mmo, 'con', 'sl1')`

By iterating this, the mmo object adds two columns per operation (log2FC and padj) for each comparison.

```
mmo <- PairwiseComp(mmo, 'con', 'sl2')
mmo <- PairwiseComp(mmo, 'con', 'sl3')
mmo <- PairwiseComp(mmo, 'con', 'px1')
mmo <- PairwiseComp(mmo, 'con', 'px2')
mmo <- PairwiseComp(mmo, 'con', 'px3')
mmo <- PairwiseComp(mmo, 'con', 'le1')
mmo <- PairwiseComp(mmo, 'con', 'le2')
mmo <- PairwiseComp(mmo, 'con', 'le3')
```
