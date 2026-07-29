# GetRichness

Sample-level richness: number of features present in each sample. A
feature is present if value \> threshold; 0 never counts as present.

## Usage

``` r
GetRichness(feature_data, metadata, threshold = 0)
```

## Arguments

  - feature\_data:
    
    Feature table with columns: id, feature, then sample columns

  - metadata:
    
    Metadata table with sample and group columns

  - threshold:
    
    Numeric; detection threshold for presence (default: 0)

## Value

data.frame with columns: sample, group, richness
