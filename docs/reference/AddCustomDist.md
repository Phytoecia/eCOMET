# Add a custom feature distance matrix to the `mmo`

Stores any user-supplied pairwise feature distance matrix in the `mmo`
so it can be used by `GetBetaDiversity()`, `GetAlphaDiversity()`, and
`HeatmapPlot()` via the `distance` argument.

## Usage

``` r
AddCustomDist(mmo, dist_matrix, name)
```

## Arguments

  - mmo:
    
    The `mmo`

  - dist\_matrix:
    
    A square numeric matrix of pairwise feature dissimilarities (see
    Details for format requirements).

  - name:
    
    A short string used to identify this distance in downstream
    functions (e.g. `"tanimoto"`, `"npc"`). The matrix is stored as
    `mmo$<name>.dissim` and referenced by passing `distance = "<name>"`
    to `GetBetaDiversity()`, `GetAlphaDiversity()`, etc. Must not
    conflict with existing names: `"dreams"`, `"cosine"`, `"m2ds"`.

## Value

The `mmo` with the new distance matrix stored in `mmo$<name>.dissim`.

## Details

**Required matrix format:**

  - A square numeric matrix with equal row and column names.

  - Row and column names must be feature IDs matching the `id` column in
    `mmo$feature_data` (e.g. `"1234.56_2.34"` or whatever identifier
    MZmine assigned).

  - Values must be **dissimilarities** in the range `[0, 1]`, where 0
    means identical and 1 means maximally different. If your source data
    are similarities (e.g. Tanimoto coefficients, cosine scores),
    convert first with `dist_matrix <- 1 - sim_matrix`.

  - The diagonal should be 0 (a feature is identical to itself).

  - The matrix should be symmetric (`dist[i,j] == dist[j,i]`).
    Non-symmetric matrices are accepted but may produce unexpected
    results in downstream ordination and Hill number calculations.

  - Features not present in `mmo$feature_data` are retained in the
    matrix but will be silently ignored during analysis. Features in
    `mmo$feature_data` that are absent from the matrix will be excluded
    from structure-aware calculations.

## Examples

``` r
if (FALSE) {
# Example: add a Tanimoto-based structural distance
# tanimoto_sim is a feature x feature similarity matrix from an external tool
tanimoto_dist <- 1 - tanimoto_sim
mmo <- AddCustomDist(mmo, dist_matrix = tanimoto_dist, name = "tanimoto")
# Now use it anywhere a distance argument is accepted:
beta <- GetBetaDiversity(mmo, method = "CSCS", distance = "tanimoto")
alpha <- GetAlphaDiversity(mmo, mode = "weighted", distance = "tanimoto")
}
```
