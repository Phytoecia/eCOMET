# GetAlphaDiversity (similarity-based)

Identical to `GetAlphaDiversity_derep()` but uses `GetSimMat()` to
retrieve similarity matrices from `.sim` slots (created by
`AddChemSim()`) and passes them to the similarity-aware inner functions
(`GetFunctionalHillNumber()`, `GetFaithPD()`).

## Usage

``` r
GetAlphaDiversity(
  mmo,
  q = 1,
  normalization = "None",
  mode = "richness",
  distance = "dreams",
  threshold = 0,
  filter_id = FALSE,
  id_list = NULL,
  filter_group = FALSE,
  group_list = NULL,
  output = c("sample_level", "group_average", "group_cumulative", "rarefied_sample"),
  group_col = "group",
  sample_col = "sample",
  pool_method = c("sum", "mean"),
  n_perm = 200,
  ci = 0.95,
  seed = NULL,
  use_mst = FALSE,
  use_fastcluster = FALSE
)
```

## Arguments

  - mmo:
    
    The `mmo` containing feature data and metadata

  - q:
    
    The Hill number order controlling abundance sensitivity (default:
    1). Only applies to `mode = "weighted"` and `mode = "unweighted"`;
    ignored for `"richness"` and `"faith"`.
    
      - `q = 0` – richness: all detected features count equally
        regardless of abundance.
    
      - `q = 1` – Shannon-type: features weighted proportionally to
        their relative abundance.
    
      - `q = 2` – Simpson-type: dominant (high-abundance) features
        weighted more strongly.

  - normalization:
    
    Abundance table to use. Options: 'None', 'Log', 'Meancentered', 'Z',
    'PA' (default: 'None'). Using `'PA'` forces presence/absence
    regardless of `mode` or `q`, effectively making every detected
    feature equally abundant before the Hill calculation.

  - mode:
    
    The diversity metric to calculate. One of 'weighted', 'unweighted',
    'faith', 'richness' (default: 'richness'). Use `q` to control
    abundance sensitivity for 'weighted' and 'unweighted'.

  - distance:
    
    Feature dissimilarity metric: 'dreams', 'm2ds', or 'cosine'
    (default: 'dreams'). Required for `mode = "weighted"` and `mode =
    "faith"`; ignored otherwise.

  - threshold:
    
    Numeric threshold used to define metabolite presence (default: 0)

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

  - output:
    
    Output mode: 'sample\_level', 'group\_average', 'group\_cumulative',
    or 'rarefied\_sample'

  - group\_col:
    
    Column in mmo$metadata that defines groups (default: 'group')

  - sample\_col:
    
    Column in mmo$metadata that defines sample IDs (default: 'sample')

  - pool\_method:
    
    How to pool abundances when combining samples: 'sum' or 'mean'
    (default: 'sum')

  - n\_perm:
    
    Integer; maximum number of permutations per rarefaction level
    (default: 200)

  - ci:
    
    Numeric; confidence level (default: 0.95)

  - seed:
    
    Optional integer seed for reproducibility (default: NULL)
