# GetFunctionalHillNumber (similarity-based)

Calculates the functional Hill number using a pairwise feature
**similarity** matrix (values 0–1, higher = more similar) from
`GetSimMat()`. Algebraically equivalent to
`GetFunctionalHillNumber_derep()` but avoids ever constructing the full
dense distance matrix, making it safe for sparse large-dataset
similarity matrices. Key identity used: \\(p'Dp = 1 - p'Sp\\) (exact
when \\(\\sum p\_i = 1\\) and \\(D = 1 - S\\)).

## Usage

``` r
GetFunctionalHillNumber(
  feature,
  metadata,
  sim_matrix,
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

  - sim\_matrix:
    
    Feature similarity matrix (dense or sparse dgCMatrix)

  - q:
    
    The order of the Hill number (default: 1)

  - threshold:
    
    Numeric; detection threshold for presence (default: 0)

  - scale\_dissim:
    
    Kept for interface compatibility; max(sim)=1 so scaling is a no-op,
    but passing FALSE still skips it

## Value

A data frame with sample, group, hill\_number, and value columns.
