# GetFaithPD (similarity-based)

Calculates Faith's phylogenetic diversity using a pairwise feature
**similarity** matrix from `GetSimMat()`. Converts to distance
internally (`1 - sim`) before clustering, which requires materializing a
full dense matrix; above 10,000 features it switches to a minimum
spanning tree, which stays sparse but uses single-linkage rather than
average-linkage.

## Usage

``` r
GetFaithPD(
  feature,
  metadata,
  sim_matrix,
  threshold = 0,
  use_mst = FALSE,
  use_fastcluster = FALSE
)
```

## Arguments

  - feature:
    
    Feature table with columns: id, feature, then sample columns

  - metadata:
    
    Metadata table with sample and group columns

  - sim\_matrix:
    
    Feature similarity matrix (dense or sparse dgCMatrix)

  - threshold:
    
    Numeric; detection threshold for presence (default: 0)

  - use\_mst:
    
    Logical; if TRUE, always use MST/single-linkage tree (works on
    sparse matrices of any size). If FALSE (default), average-linkage is
    used for n \<= 10,000 and MST is used automatically (with a warning)
    for larger matrices.

  - use\_fastcluster:
    
    Logical; if TRUE, build the average-linkage tree with
    fastcluster::hclust() instead of stats::hclust(). Both are valid
    UPGMA implementations, but they break ties between equal distances
    differently. Chemical similarity scores are usually reported to
    three decimals, so ties are common and the two engines can return
    different trees. The default FALSE keeps results identical to
    earlier releases and to machines that do not have fastcluster
    installed.

## Value

A data frame with sample, group, PD, SR, and value columns.
