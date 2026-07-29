# Annotate `mmo$feature_info` with MS2 presence and MS2 block counts from an MGF

Scan an `.mgf` file and summarize MS/MS availability for each feature in
`mmo$feature_info`. The function adds two columns:

  - `ms2`: `TRUE` if the MGF contains at least one `MSLEVEL=2` block for
    that `id`; otherwise `FALSE`.

  - `count_ms2`: number of `MSLEVEL=2` blocks for that `id`.

## Usage

``` r
annotate_feature_info_ms2_from_mgf(
  mmo,
  mgf_path,
  chunk_lines = 100000L,
  overwrite = FALSE,
  verbose = TRUE
)
```

## Arguments

  - mmo:
    
    An eCOMET `mmo` containing a required `feature_info` table with an
    `id` column (`mmo$feature_info$id`).

  - mgf\_path:
    
    Character. Path to the input `.mgf` file.

  - chunk\_lines:
    
    Integer. Number of lines read per iteration. Larger values are
    typically faster but use more memory. Default is `100000L`.

  - overwrite:
    
    Logical. If `FALSE` (default) and `ms2` and/or `count_ms2` already
    exist in `mmo$feature_info`, the function errors. Set `overwrite =
    TRUE` to replace existing columns.

  - verbose:
    
    Logical. If `TRUE` (default), prints a brief summary of how many MS2
    blocks were found and how many features have MS2.

## Value

The updated `mmo` with `mmo$feature_info$ms2` and
`mmo$feature_info$count_ms2` added (or overwritten if `overwrite =
TRUE`).

## Details

This is useful for quickly identifying which features have MS/MS spectra
available (and how many replicate MS2 spectra exist) before downstream
annotation, networking, or library-building steps.

## Examples

``` r
if (FALSE) { # \dontrun{
# Add ms2 + count_ms2 columns to mmo$feature_info
mmo <- annotate_feature_info_ms2_from_mgf(mmo,
       "spectra.mgf")

# Overwrite existing columns if you re-run on a different MGF
mmo <- annotate_feature_info_ms2_from_mgf(mmo,
       "spectra_new.mgf", overwrite = TRUE)
} # }
```
