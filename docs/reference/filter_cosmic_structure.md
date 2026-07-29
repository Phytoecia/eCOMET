# Filter SIRIUS structure (CSI:FingerID) annotations by COSMIC confidence score

Applies a minimum COSMIC confidence cutoff to SIRIUS structure
predictions inside an ecomet `mmo`. The function reads a chosen
annotation table from `mmo[[input]]` (default `"sirius_annot"`), flags
structure annotations below `threshold` by setting selected
structure-identification fields to `NA`, and stores the result as a new
element on `mmo` named `"sirius_annot_filtered_<suffix>"`.

## Usage

``` r
filter_cosmic_structure(
  mmo,
  input = "sirius_annot",
  cosmic_mode = c("exact", "approx"),
  threshold = 0.5,
  fields = "auto",
  na_cosmic = TRUE,
  suffix = NULL,
  overwrite = FALSE,
  verbose = TRUE
)
```

## Arguments

  - mmo:
    
    An eCOMET `mmo` containing `mmo[[input]]` (a data.frame).

  - input:
    
    Character. Name of the element on `mmo` to filter. Defaults to
    `"sirius_annot"`.

  - cosmic\_mode:
    
    Which COSMIC column to use. One of `"exact"` or `"approx"`.

  - threshold:
    
    Numeric. Keep structures with COSMIC \>= threshold.

  - fields:
    
    Character vector of columns to NA-out when COSMIC \< threshold. If
    `"auto"` (default), uses a reasonable default set if present in the
    table.

  - na\_cosmic:
    
    Logical. If `TRUE` (default), also set the COSMIC value to `NA` when
    the structure is filtered out.

  - suffix:
    
    Optional character string appended to the created element name.

  - overwrite:
    
    Logical. If `FALSE` (default) and the target element already exists,
    error.

  - verbose:
    
    Logical. If `TRUE`, prints a concise summary.

## Value

The updated `mmo`, with a new element
`mmo[[paste0("sirius_annot_filtered_", suffix)]]`.

## Details

Rows are never dropped. "Removing" a structure means setting selected
fields (e.g., SMILES/InChI/name/InChIkey2D) to `NA`.

## Examples

``` r
if (FALSE) { # \dontrun{
mmo <- filter_cosmic_structure(mmo, input = "sirius_annot",
                              cosmic_mode = "approx", threshold = 0.5,
                              suffix = "COSMIC_exact_0.5", verbose = TRUE)
} # }
```
