# Add a custom feature similarity matrix to the `mmo`

Stores any user-supplied pairwise feature similarity matrix in the `mmo`
so it can be used by `GetBetaDiversity()`, `GetAlphaDiversity()`, and
`HeatmapPlot()` via the `distance` argument. Accepts either a dense
numeric matrix or a sparse `dgCMatrix` (from the Matrix package).

## Usage

``` r
AddCustomSim(mmo, sim_matrix, name)
```

## Arguments

  - mmo:
    
    The `mmo`

  - sim\_matrix:
    
    A square numeric matrix or sparse Matrix of pairwise feature
    similarities (see Details for format requirements).

  - name:
    
    A short string used to identify this matrix in downstream functions
    (e.g. `"tanimoto"`, `"npc"`). Stored as `mmo$<name>.sim` and
    referenced via `distance = "<name>"`. Must not conflict with
    built-in names: `"dreams"`, `"cosine"`, `"m2ds"`.

## Value

The `mmo` with the similarity matrix stored in `mmo$<name>.sim`.

## Details

**Required matrix format:**

  - A square numeric matrix (or sparse `dgCMatrix`) with equal row and
    column names.

  - Row and column names must be feature IDs matching the `id` column in
    `mmo$feature_data`.

  - Values must be **similarities** in the range `[0, 1]`, where 1 means
    identical and 0 means maximally different.

  - The diagonal should be 1 (a feature is identical to itself).

  - The matrix should be symmetric (`sim[i,j] == sim[j,i]`).

  - Features not present in `mmo$feature_data` are retained but silently
    ignored during analysis.

## Examples

``` r
if (FALSE) {
# tanimoto_sim is a feature x feature similarity matrix from an external tool
mmo <- AddCustomSim(mmo, sim_matrix = tanimoto_sim, name = "tanimoto")
# Now use it anywhere a distance argument is accepted:
beta <- GetBetaDiversity(mmo, method = "CSCS", distance = "tanimoto")
alpha <- GetAlphaDiversity(mmo, mode = "weighted", distance = "tanimoto")
}
```
