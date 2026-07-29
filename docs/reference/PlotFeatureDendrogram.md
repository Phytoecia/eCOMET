# PlotFeatureDendrogram

Plot a feature dendrogram produced by `FeatureDendrogram()` using
ggtree, returning a `ggplot`/`ggtree` object that can be further
customised with additional `ggplot2` or `ggtree` layers.

## Usage

``` r
PlotFeatureDendrogram(
  tree,
  mmo,
  color_by = "NPC#pathway",
  layout = "circular",
  color_branches = FALSE,
  highlight_groups = TRUE,
  palette = "Dark 3",
  na_color = "grey70",
  tip_size = 1,
  branch_size = 0.5,
  show_tip_labels = FALSE,
  open_angle = 10
)
```

## Arguments

  - tree:
    
    Output list from `FeatureDendrogram()`.

  - mmo:
    
    `mmo`. Used for `mmo$sirius_annot` when `color_by` is not `NULL`.

  - color\_by:
    
    Column name in `mmo$sirius_annot` to use for tip point colors
    (default: `"NPC#pathway"`). Set to `NULL` to skip tip coloring.

  - layout:
    
    Tree layout passed to `ggtree::ggtree()`. One of `"circular"`
    (default), `"rectangular"`, `"fan"`, `"slanted"`, `"radial"`, or any
    other layout supported by ggtree.

  - color\_branches:
    
    Logical; whether to color the tree branches themselves by the
    annotation in `color_by`, in addition to the tip points (default:
    `FALSE`). Branch color is determined by a majority-vote of the tip
    annotations descending from each node, so mixed-class clades are
    colored by whichever class is most abundant. Unannotated or mixed
    branches fall back to `na_color`.

  - highlight\_groups:
    
    Logical; whether to highlight ion identity / correlation groups
    using `ggtree::geom_hilight()` (default: `TRUE`). Only has an effect
    when `tree$tip_map` is not `NULL`.

  - palette:
    
    Qualitative palette name for `colorspace::qualitative_hcl` (default:
    `"Dark 3"`).

  - na\_color:
    
    Color for tips/branches with no annotation (default: `"grey70"`).

  - tip\_size:
    
    Point size for tip points (default: `1`).

  - branch\_size:
    
    Line width for tree branches (default: `0.5`). Increase for small
    trees; decrease for large ones.

  - show\_tip\_labels:
    
    Logical; draw feature ID text at tips (default: `FALSE`). Only
    useful for small trees.

  - open\_angle:
    
    For `layout = "fan"`, the opening angle in degrees (default: `10`).

## Value

A `ggtree`/`ggplot` object. Add further layers with `+`.

## Details

**Design philosophy:** The function returns a `ggtree` object so you can
keep layering modifications after the call:

``` 
  p <- PlotFeatureDendrogram(tree, mmo)
  p + ggtree::geom_tiplab(size = 1.5) +
      ggplot2::theme(legend.position = "bottom")
```

**Tip point coloring (`color_by`):** Because most metabolomics trees
have hundreds or thousands of tips, text labels are rarely readable.
Instead, tiny colored points are drawn at each tip using
`ggtree::geom_tippoint()` to indicate compound class. Use
`show_tip_labels = TRUE` to additionally draw feature ID text (not
recommended for large trees).

**IIN / correlation group highlighting (`highlight_groups`):** When
`tree$tip_map` is present, each ion identity / correlation group is
highlighted using `ggtree::geom_hilight()`, which shades the MRCA clade
of all members of that group.

## Examples

``` r
if (FALSE) {
tree <- FeatureDendrogram(mmo, distance = "dreams")

# Basic circular tree coloured by NPC pathway (tip points only)
PlotFeatureDendrogram(tree, mmo)

# Color branches too, thicker lines for a small tree
PlotFeatureDendrogram(tree, mmo, color_branches = TRUE, branch_size = 1.2)

# Rectangular layout, no coloring
PlotFeatureDendrogram(tree, mmo, layout = "rectangular", color_by = NULL)

# Returned object is a ggplot -- keep layering
p <- PlotFeatureDendrogram(tree, mmo)
p + ggtree::geom_tiplab(size = 1, align = TRUE)
}
```
