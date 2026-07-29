# ExportITOL

Export a feature dendrogram and companion annotation files for
visualisation in iTOL (Interactive Tree of Life). Produces three files:

1.  A Newick tree (`<outprefix>.nwk`)

2.  A colour-strip annotation coloured by NPC pathway (or any
    `mmo$sirius_annot` column) (`<outprefix>_colorstrip.txt`)

3.  A bar-chart annotation showing per-feature prevalence across samples
    (`<outprefix>_barplot.txt`)

## Usage

``` r
ExportITOL(
  tree,
  mmo,
  outprefix = "itol_export",
  color_by = "NPC#pathway",
  palette = "Dark 3",
  na_color = "#CCCCCC"
)
```

## Arguments

  - tree:
    
    Output list from `FeatureDendrogram()`.

  - mmo:
    
    `mmo`. Must contain `mmo$sirius_annot` (for colour strip) and
    `mmo$feature_presence` or `mmo$feature_data` (for prevalence bar
    chart).

  - outprefix:
    
    File path prefix for output files (default: `"itol_export"`).
    Directories in the path must exist.

  - color\_by:
    
    Column name in `mmo$sirius_annot` used for the colour strip
    (default: `"NPC#pathway"`).

  - palette:
    
    Qualitative palette name for `colorspace::qualitative_hcl` (default:
    `"Dark 3"`).

  - na\_color:
    
    Hex colour for features with no annotation (default: `"#CCCCCC"`).

## Value

Invisibly returns a named character vector of written file paths.

## Details

**How to load in iTOL:** Upload the Newick file at
<https://itol.embl.de/upload.cgi>, then drag and drop the two annotation
files onto the displayed tree. Select "Circular" layout and enable
"Display" for each dataset in the Datasets panel.

**Prevalence calculation:** Feature prevalence is the proportion of
samples in which a feature is detected (abundance \> 0). Features absent
from the PA table are assigned prevalence 0. This value drives the
bar-chart length in iTOL.

**Colour-strip annotation:** Each tip is assigned the colour of its NPC
pathway (or other chosen column) from `mmo$sirius_annot`. Features with
no annotation are coloured `na_color`. The legend is written into the
file header so iTOL renders it automatically.

## Examples

``` r
if (FALSE) {
tree <- FeatureDendrogram(mmo, distance = "dreams")
ExportITOL(tree, mmo, outprefix = "output/itol_dreams")
# Upload output/itol_dreams.nwk to iTOL, then drag and drop the two .txt files
}
```
