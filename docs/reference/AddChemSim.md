# Add chemical similarity matrices to the `mmo`

This function reads cosine, DREAMS, and MS2DeepScore molecular
networking outputs from MZmine and stores the raw pairwise
**similarity** matrices (values 0–1, higher = more similar) in the
`mmo`. For datasets with more than 10,000 features a sparse `Matrix` is
used to avoid allocating a full dense matrix; unobserved feature pairs
carry implicit similarity 0 (equivalent to dissimilarity 1, maximally
different).

## Usage

``` r
AddChemSim(mmo, cos_dir = NULL, dreams_dir = NULL, m2ds_dir = NULL)
```

## Arguments

  - mmo:
    
    The `mmo`

  - cos\_dir:
    
    Path to the cosine similarity CSV file from MZMine (molecular
    networking)

  - dreams\_dir:
    
    Path to the DREAMS similarity CSV file from MZMine (molecular
    networking)

  - m2ds\_dir:
    
    Path to the MS2DeepScore similarity CSV file from MZMine (molecular
    networking)

## Value

The `mmo` with similarity matrices stored in `mmo$cos.sim`,
`mmo$dreams.sim`, and/or `mmo$m2ds.sim`. Large datasets (n \> 10,000)
use sparse `dgCMatrix` format; smaller datasets use a dense numeric
matrix.

## Note

Tree-based downstream functions (`GetFaithPD`, `FeatureDendrogram`,
`GetBetaDiversity` with `method = "Gen.Uni"`) require a full n x n
distance matrix and will fail at large scale regardless of sparse
storage. Use `filter_mmo()` to reduce feature count before calling those
functions.

## Examples

``` r
if (FALSE) {
mmo <- AddChemSim(mmo,
 cos_dir = "path/to/cosine_similarity.csv",
 dreams_dir = "path/to/dreams_similarity.csv",
 m2ds_dir = "path/to/ms2deepscore_similarity.csv"
)
}
```
