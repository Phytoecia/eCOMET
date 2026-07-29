# Adding annotation from SIRIUS to the `mmo`

This function reads SIRIUS structure identification and formula summary
files, and adds the annotations to the `mmo`.

## Usage

``` r
AddSiriusAnnot(
  mmo,
  canopus_structuredir,
  canopus_formuladir,
  filter_annot = FALSE,
  filter_threshold = 0.5
)
```

## Arguments

  - mmo:
    
    The `mmo`

  - canopus\_structuredir:
    
    Path to the SIRIUS structure\_identification.tsv file

  - canopus\_formuladir:
    
    Path to the SIRIUS canopus\_formula\_summary.tsv file

  - filter\_annot:
    
    Logical. If TRUE, filter the annotations by probability threshold in
    CANOPUS.

  - filter\_threshold:
    
    Numeric between 0 and 1. The probability threshold for filtering
    annotations.

## Value

The `mmo` with SIRIUS annotations added

## Examples

``` r
if (FALSE) {
mmo <- AddSiriusAnnot(mmo,
 canopus_structuredir = "path/to/structure_identification.tsv",
 canopus_formuladir = "path/to/canopus_formula_summary.tsv"
)
}
```
