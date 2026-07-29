# FeatureDendrogram

Build a feature-to-feature dendrogram from a dissimilarity matrix stored
in the `mmo`. Optionally collapses within-group distances for ion
identity or correlation groups so that adducts/isotopes of the same
compound cluster together before chemical distance drives the
higher-level topology.

## Usage

``` r
FeatureDendrogram_derep(
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

  - `dist_used` – the (possibly modified) distance matrix used

  - `tip_map` – data.frame of feature ID -\> group assignment (`NULL`
    when `ion_identity = "none"`)

## Details

**Ion identity options:**

  - `"none"` – use the raw dissimilarity matrix as-is.

  - `"correlation"` – features sharing the same MZmine correlation group
    (typically co-eluting features likely derived from the same
    compound) are assigned `within_group_dist` before clustering.

  - `"ion_identity_network"` – features linked by MZmine ion identity
    networking (confirmed adducts / isotopes) are assigned
    `within_group_dist` before clustering.

In both grouped modes the constraint is applied symmetrically and the
diagonal is forced to 0 after the assignment.

## Examples

``` r
if (FALSE) {
tree <- FeatureDendrogram(mmo, distance = "dreams")
tree_iin <- FeatureDendrogram(mmo, distance = "dreams",
                              ion_identity = "ion_identity_network")
}
```
