# GetBetaDiversity

Calculate beta diversity (sample-to-sample dissimilarity) for a `mmo`.
Four methods are supported, differing in how they handle feature
abundance and structural relationships among compounds:

## Usage

``` r
GetBetaDiversity_derep(
  mmo,
  method = "Gen.Uni",
  normalization = "None",
  distance = NULL,
  filter_id = FALSE,
  id_list = NULL,
  filter_group = FALSE,
  group_list = NULL,
  scale_dissim = TRUE
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

## Value

For 'bray', 'jaccard', 'CSCS': a symmetric sample-by-sample distance
matrix (matrix). For 'Gen.Uni': a named list with three distance
matrices named `d_0` (presence/absence), `d_0.5` (balanced), and `d_1`
(fully abundance-weighted).

## Details

  - 'bray' : Bray-Curtis dissimilarity. Feature-overlap only; no
    compound distance matrix required. Abundance sensitivity is
    controlled via `normalization`: use `'None'` for raw abundances,
    `'Log'` to compress dynamic range, or `'PA'` to reduce to
    presence/absence.

  - 'jaccard' : Jaccard dissimilarity. Presence/absence only by
    definition. The `normalization` argument is ignored; PA data are
    always used internally. A warning is issued if a non-PA
    normalization is supplied.

  - 'CSCS' : Chemical Structural and Compositional Similarity.
    Incorporates pairwise compound distances so that samples sharing
    structurally related (but not identical) features can still be
    considered similar. Requires `distance`. Abundance sensitivity is
    controlled via `normalization` (same options as bray).

  - 'Gen.Uni' : Generalized UniFrac. Partitions beta diversity across
    branches of a compound dendrogram. Requires `distance`. Returns a
    named list of three distance matrices, one per abundance-weighting
    level (`d_0`, `d_0.5`, `d_1`). Choose the one appropriate for your
    analysis (see Details).

**Choosing an abundance-weighting level for Gen.Uni:**

  - `d_0` – presence/absence only; equivalent to unweighted UniFrac. Use
    when detection (not intensity) is the signal of interest.

  - `d_0.5` – balanced weighting. Moderates the influence of highly
    abundant features without ignoring abundance.

  - `d_1` – fully abundance-weighted. Dominant features drive the
    distance. Use when peak intensity is a reliable biological signal.

Access a specific matrix with e.g. `result[["d_0.5"]]`.

## Examples

``` r
if (FALSE) {
beta_diversity <- GetBetaDiversity(mmo, method = 'Gen.Uni',
 normalization = 'None', distance = 'dreams', filter_id = FALSE)
# Access a specific weighting level:
d_balanced <- beta_diversity[["d_0.5"]]
beta_diversity <- GetBetaDiversity(mmo, method = 'bray',
 normalization = 'Z', filter_id = TRUE, id_list = Glucosinolates,
 filter_group = TRUE, group_list = c('Control', 'Treatment1'))
}
```
