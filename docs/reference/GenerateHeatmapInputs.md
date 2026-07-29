# GenerateHeatmapInputs (similarity-based)

Identical to `GenerateHeatmapInputs_derep()` but retrieves feature
similarity matrices via `GetSimMat()` (.sim slots) and converts to
distance (`1 - S`) when building the `dist` object for heatmap row
clustering.

## Usage

``` r
GenerateHeatmapInputs(
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

  - dist\_matrix: A distance matrix derived from the similarity matrix
    as `1 - S`, or `NULL` when `distance` is not supplied

  - row\_label: A vector of row labels for custom-annotated features
    (see `AddCustomAnnot()`). Feature IDs are used when no custom
    annotation is available.

  - heatmap\_data: A data frame of the heatmap values with feature IDs
