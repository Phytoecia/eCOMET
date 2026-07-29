# Filter an MGF file to keep only spectra for features present in `mmo$feature_data$id`

Create a new `.mgf` file that contains only spectra (ION blocks) whose
`FEATURE_ID` occurs in `mmo$feature_data$id`. This is useful for keeping
your spectral library in sync with the features currently stored in a
`mmo` (e.g., after subsetting, filtering, or rebuilding a `mmo`).

## Usage

``` r
filter_mgf_to_mmo(
  mmo,
  mgf_path,
  output_path = NULL,
  chunk_lines = 100000L,
  verbose = TRUE
)
```

## Arguments

  - mmo:
    
    An eCOMET `mmo` containing a required `feature_data` table with an
    `id` column (`mmo$feature_data$id`).

  - mgf\_path:
    
    Character. Path to the input `.mgf` file.

  - output\_path:
    
    Character or NULL. Path to write the filtered `.mgf`. If `NULL`
    (default), the output name is derived from `mgf_path` by appending
    `_filtered` before the `.mgf` extension (or adding `_filtered.mgf`
    if no extension is present).

  - chunk\_lines:
    
    Integer. Number of lines read per iteration. Larger values are
    typically faster but use more memory. Default is `100000L`.

  - verbose:
    
    Logical. If `TRUE` (default), prints a short progress summary and
    final counts.

## Value

Invisibly returns a list with:

  - `output_path`: path to the filtered MGF

  - `blocks_total`: total `BEGIN IONS` blocks encountered

  - `blocks_kept`: number of blocks written to `output_path`

  - `lines_read`: total lines read from `mgf_path`

## Details

Each spectrum in an MGF is represented by a `BEGIN IONS ... END IONS`
block. This function keeps or discards entire blocks based on the
integer value in the `FEATURE_ID=` header line. All blocks (MS1 and MS2)
are kept for a retained feature, including multiple MS2 blocks if
present.

## Examples

``` r
if (FALSE) { # \dontrun{
# Write "<input>_filtered.mgf" containing only FEATURE_IDs present in mmo$feature_data$id
filter_mgf_to_mmo(mmo, "spectra.mgf")

# Write to a custom file name
filter_mgf_to_mmo(mmo, "spectra.mgf", output_path = "spectra_mmo_only.mgf")
} # }
```
