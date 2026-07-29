# GetBetaDiversity (similarity-based)

Identical to `GetBetaDiversity_derep()` but retrieves pairwise
**similarity** matrices via `GetSimMat()` (.sim slots). Key differences
from the \_derep version:

  - CSCS: uses the similarity matrix directly as CSS (no `1 - D`
    conversion needed), because `CSS = S` when similarities are stored.

  - Gen.Uni: converts `S -> D` internally (`1 - S`) only for the
    `hclust()` call; guarded by a size check (n \> 10,000 errors).

  - bray / jaccard: unchanged (no distance matrix used).

## Usage

``` r
GetBetaDiversity(
  mmo,
  method = "Gen.Uni",
  normalization = "None",
  distance = NULL,
  filter_id = FALSE,
  id_list = NULL,
  filter_group = FALSE,
  group_list = NULL,
  scale_dissim = TRUE,
  use_mst = FALSE,
  use_fastcluster = FALSE
)
```

## Arguments

  - mmo:
    
    The `mmo` containing feature data and metadata

  - method:
    
    Beta diversity method: 'Gen.Uni', 'bray', 'jaccard', or 'CSCS'
    (default: 'Gen.Uni')

  - normalization:
    
    Abundance table to use. Options: 'None', 'Log', 'Meancentered', 'Z',
    'PA' (default: 'None'). Ignored for 'jaccard', which always uses PA.
    For 'bray' and 'CSCS', this is the primary lever for controlling
    abundance sensitivity.

  - distance:
    
    Feature dissimilarity metric: 'dreams', 'm2ds', or 'cosine'.
    Required for 'Gen.Uni' and 'CSCS'; ignored for 'bray' and 'jaccard'.

  - filter\_id:
    
    A boolean indicating whether to filter the feature data by a
    specific list (default: FALSE)

  - id\_list:
    
    A list of feature names to filter the feature data by, if filter\_id
    is TRUE (default: NULL)

  - filter\_group:
    
    A boolean indicating whether to filter the feature data by a
    specific group list (default: FALSE)

  - group\_list:
    
    A list of groups to filter the feature data by, if filter\_group is
    TRUE (default: NULL)

  - scale\_dissim:
    
    Boolean; whether to scale the feature distance matrix to between 0,1
    (default: TRUE)

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
