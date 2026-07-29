# Generate input files to be used for pheatmap from the `mmo`

This function generates heatmap inputs from the `mmo`, including fold
change or mean values, distance matrix, and row labels for
custom-annotated features.

## Usage

``` r
GenerateHeatmapInputs_derep(
  mmo,
  filter_id = FALSE,
  id_list = NULL,
  filter_group = FALSE,
  group_list = NULL,
  summarize = "mean",
  control_group = "ctrl",
  normalization = "None",
  distance = NULL
)
```

## Arguments

  - mmo:
    
    The `mmo` with sirius annotation and normalized data

  - filter\_id:
    
    Boolean to filter features by id\_list (default: FALSE)

  - id\_list:
    
    A vector of feature names to filter (default: NULL)

  - filter\_group:
    
    Boolean to filter groups by group\_list (default: FALSE)

  - group\_list:
    
    A vector of group names to filter (default: NULL)

  - summarize:
    
    The summarization method to use. Options are 'fold\_change' or
    'mean' (default: 'mean')

  - control\_group:
    
    The group to use as control for fold change calculation (default:
    'ctrl')

  - normalization:
    
    The normalization method to use. Options are 'None', 'Log',
    'Meancentered', or 'Z'

  - distance:
    
    The distance metric to use. Options are 'dreams', 'cosine', or
    'm2ds' (default: 'dreams')

## Value

A list containing the following elements:

  - FC\_matrix: A matrix of fold change or mean values

  - dist\_matrix: A distance matrix based on the specified distance
    metric

  - row\_label: A vector of row labels for custom-annotated features
    (See AddCustomAnnot()). If no custom annotation is available,
    feature IDs are used.

  - heatmap\_data: A data frame containing the heatmap data with feature
    IDs and values

## Examples

``` r
if (FALSE) {
# Generate heatmap inputs to visualize fold change values with log normalization and dreams distance
heatmap_inputs <- GenerateHeatmapInputs(
 mmo, summarize = 'fold_change', control_group = 'Control',
 normalization = 'None', distance = 'dreams'
)
# Generate heatmap inputs to visualize mean values
heatmap_inputs <- GenerateHeatmapInputs(
 mmo, summarize = 'mean', normalization = 'None', distance = 'dreams'
)
# The resulting list contains FC_matrix, dist_matrix, row_label, and heatmap_data
# A heatmap can be generated using pheatmap
# 'clustering_distance_rows' option make the dendrogram follows chemical distances of features.
#  -Delete this option to visualize the heatmap following cannonical clustering
pheatmap(mat = heatmap_inputs$FC_matrix,
    cluster_rows = TRUE, #do not change
    clustering_distance_rows = heatmap_inputs$dist_matrix,
    cluster_cols = TRUE,
    clustering_method = "average", #UPGMA
    show_rownames = TRUE,
    show_colnames = TRUE,
    cellwidth = 25,
    cellheight = 0.05,
    treeheight_row = 100,
    fontsize_row = 3,
    fontsize_col = 15,
    scale = 'none',
    annotation_names_row = TRUE,
    labels_row = heatmap_inputs$row_label,
    )
}
```
