# FeatureDendrogram (similarity-based)

Identical to `FeatureDendrogram_derep()` but retrieves the feature
**similarity** matrix via `GetSimMat()` (.sim slots) and converts to
distance (`1 - S`) only inside `hclust()`. Adds a size guard: stops with
a clear message if n \> 10,000 features, because `as.dist()` on an n x n
matrix is infeasible at that scale. Ion identity constraints are
expressed as high similarity (`1 - within_group_dist`) rather than low
distance.

## Usage

``` r
FeatureDendrogram(
  mmo,
  distance = "dreams",
  features = NULL,
  method = "average",
  ion_identity = c("none", "correlation", "ion_identity_network"),
  corr_col = "feature_group",
  iin_col = "ion_identities:iin_id",
  within_group_dist = 0.01,
  save_newick = FALSE,
  outprefix = "feature_tree"
)
```

## Arguments

  - mmo:
    
    `mmo`. Must contain the requested dissimilarity matrix (added via
    `AddChemDist()` or `AddCustomDist()`) and, for ion identity modes,
    `mmo$feature_info`.

  - distance:
    
    Name of the dissimilarity matrix to use (default: `"dreams"`).
    Passed to `GetDistanceMat()`; supports built-in and custom names.

  - features:
    
    Optional character vector of feature IDs to include. `NULL`
    (default) uses all features in the distance matrix.

  - method:
    
    hclust linkage method (default: `"average"`). Other sensible
    choices: `"complete"`, `"ward.D2"`.

  - ion\_identity:
    
    One of `"none"`, `"correlation"`, `"ion_identity_network"` (default:
    `"none"`).

  - corr\_col:
    
    Column in `mmo$feature_info` containing MZmine correlation group IDs
    (default: `"feature_group"`). Used when `ion_identity =
    "correlation"`.

  - iin\_col:
    
    Column in `mmo$feature_info` containing MZmine IIN IDs (default:
    `"ion_identities:iin_id"`). Used when `ion_identity =
    "ion_identity_network"`.

  - within\_group\_dist:
    
    Distance assigned to pairs within the same ion identity group
    (default: `0.01`). A small positive value rather than 0 avoids
    degenerate zero-height clusters in the tree while still pulling
    grouped features together. Must be in `[0, 1]`.

  - save\_newick:
    
    Logical; if `TRUE` writes a Newick file (`<outprefix>.nwk`) of the
    tree (default: `FALSE`).

  - outprefix:
    
    File prefix used when `save_newick = TRUE` (default:
    `"feature_tree"`).

## Value

A named list with:

  - `hclust` – the `hclust` object

  - `dendrogram` – the `dendrogram` object

  - `phylo` – the `phylo` object (for ape/iTOL)

  - `dist_used` – the distance matrix used for clustering (`1 -
    sim_used`)

  - `sim_used` – the (possibly modified) similarity matrix

  - `tip_map` – data.frame of feature ID -\> group assignment (`NULL`
    when `ion_identity = "none"`)
