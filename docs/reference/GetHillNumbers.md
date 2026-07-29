# GetHillNumbers

This function calculates the Hill numbers for a given `mmo`,
normalization method, and order of the Hill number without considering
feature dissimilarity.

## Usage

``` r
GetHillNumbers(feature, metadata, q = 0, threshold = 0)
```

## Arguments

  - feature:
    
    Feature table with columns: id, feature, then sample columns

  - metadata:
    
    Metadata table with sample and group columns

  - q:
    
    The order of the Hill number to calculate (default: 0)

  - threshold:
    
    Numeric; detection threshold for presence (default: 0)

## Value

A data frame containing the Hill number for each group in the metadata,
with columns for group and hill number.
