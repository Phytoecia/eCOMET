# GetFunctionalHillNumber

This function calculates the functional Hill number for a given `mmo`,
normalization method, and distance metric. See
https://nph.onlinelibrary.wiley.com/doi/full/10.1111/nph.18685 for
details of the functional Hill number calculation.

## Usage

``` r
GetFunctionalHillNumber_derep(
  feature,
  metadata,
  distance_matrix,
  q = 1,
  threshold = 0,
  scale_dissim = TRUE
)
```

## Arguments

  - feature:
    
    Feature table with columns: id, feature, then sample columns

  - metadata:
    
    Metadata table with sample and group columns

  - distance\_matrix:
    
    Feature distance matrix

  - q:
    
    The order of the Hill number to calculate (default: 1). Larger q
    values give more weight to evenness portion of the hill number over
    richness. q = 0 treats abundance matrix as a presence absence matrix

  - threshold:
    
    Numeric; detection threshold for presence (default: 0)

  - scale\_dissim:
    
    Boolean; whether to scale the distance matrix to be between 0 and 1
    (default: TRUE)

## Value

A data frame containing the functional Hill number for each group in the
metadata, with columns for group and hill number.
