# Retrieve a feature similarity matrix from the `mmo`

Looks up a similarity matrix stored in the `mmo` by name. Works with the
three built-in matrices added by `AddChemSim()` as well as any custom
matrix added via `AddCustomSim()`. Returns the matrix as-is (dense
numeric or sparse dgCMatrix); downstream functions are responsible for
converting to distance when needed.

## Usage

``` r
GetSimMat(mmo, distance = "dreams")
```

## Arguments

  - mmo:
    
    The `mmo`

  - distance:
    
    Name of the similarity matrix to retrieve. Built-in options are
    `'dreams'`, `'cosine'`, and `'m2ds'`. Any name passed to
    `AddCustomSim(mmo, name = ...)` is also valid.

## Value

The similarity matrix (numeric matrix or sparse dgCMatrix with feature
ID row/col names).

## Examples

``` r
if (FALSE) {
sim_matrix <- GetSimMat(mmo, distance = 'dreams')
sim_matrix <- GetSimMat(mmo, distance = 'tanimoto')  # custom
}
```
