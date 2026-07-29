# Filter CANOPUS / SIRIUS annotations in an eCOMET `mmo` by probability threshold

Applies a minimum-probability cutoff to selected CANOPUS (NPClassifier /
ClassyFire) annotation levels inside an eCOMET `mmo`. The function reads
a chosen annotation table from `mmo[[input]]` (default
`"sirius_annot"`), flags annotations below `threshold` by setting them
to `NA`, and stores the result as a new element on `mmo` named
`"sirius_annot_filtered_<suffix>"`.

## Usage

``` r
filter_canopus_annotations(
  mmo,
  input = "sirius_annot",
  pathway_level = "NPC#pathway",
  threshold = 0.9,
  na_prob = TRUE,
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

  - pathway\_level:
    
    Character vector of one or more annotation header(s) to filter.
    Valid options include:
    
      - `"All"`
    
      - `"All_NPC"`
    
      - `"All_ClassyFire"`
    
      - `"NPC#pathway"`
    
      - `"NPC#superclass"`
    
      - `"NPC#class"`
    
      - `"ClassyFire#superclass"`
    
      - `"ClassyFire#class"`
    
      - `"ClassyFire#subclass"`
    
      - `"ClassyFire#level 5"`
    
      - `"ClassyFire#most specific class"`
    
    Special values:
    
      - `"All"` expands to all NPC + ClassyFire levels listed above.
    
      - `"All_NPC"` expands to NPC levels only.
    
      - `"All_ClassyFire"` expands to ClassyFire levels only.

  - threshold:
    
    Decimal between 0 and 1. Annotations with probability \< threshold
    are flagged to `NA`.

  - na\_prob:
    
    Logical. If `TRUE` (default), also set the corresponding probability
    value to `NA`.

  - suffix:
    
    Optional character string appended to the created element name:
    `"sirius_annot_filtered_<suffix>"`. If `NULL` (default), a suffix is
    auto-generated.

  - overwrite:
    
    Logical. If `FALSE` (default) and the target element already exists
    on `mmo`, the function errors to avoid accidental overwrites.

  - verbose:
    
    Logical. If `TRUE`, prints a concise summary including counts of
    non-missing annotations before and after filtering.

## Value

The updated `mmo`, with a new element
`mmo[[paste0("sirius_annot_filtered_", suffix)]]`.

## Details

Rows are never dropped. "Removing" an annotation means setting the
annotation value (and optionally its probability column) to `NA`.

For each annotation column (e.g., `"NPC#pathway"`), the function looks
for an associated probability column using common SIRIUS export naming:

  - `"<header> Probability"`

  - `"<header> probability"`

  - `"<header>Probability"`

## Examples

``` r
if (FALSE) { # \dontrun{
mmo <- filter_canopus_annotations(mmo, input = "sirius_annot",
                                 pathway_level = "NPC#pathway", threshold = 0.9,
                                 suffix = "NPC_pathway_0.9", verbose = TRUE)

mmo <- filter_canopus_annotations(mmo, input = "sirius_annot_filtered_NPC_pathway_0.9",
                                 pathway_level = "All_NPC", threshold = 0.95,
                                 suffix = "All_NPC_0.95", verbose = TRUE)
} # }
```
