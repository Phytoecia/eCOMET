# GetFaithPD

This function calculates the Faith's phylogenetic diversity for a given
`mmo` and distance metric, to calculate chemically-informed richness

## Usage

``` r
GetFaithPD_derep(feature, metadata, distance_matrix, threshold = 0)
```

## Arguments

  - feature:
    
    Feature table with columns: id, feature, then sample columns

  - metadata:
    
    Metadata table with sample and group columns

  - distance\_matrix:
    
    Feature distance matrix

  - threshold:
    
    Numeric; detection threshold for presence (default: 0)

## Value

A data frame containing the Faith's phylogenetic diversity for each
group in the metadata, with columns for group and PD.
