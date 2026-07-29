# ecomet (development version)

## Similarity-based storage, review follow-ups

* `ScreenFeaturePhenotypeCorrelation()`, `GetPerformanceFeatureCorrelation()`,
  and `GetPerformanceFeatureRegression()` no longer report a correlation of
  zero when a value is missing. Samples with no phenotype value are dropped
  before the vectorised covariance step, and features carrying a missing value
  are computed with `cor.test()` so that pairwise deletion still applies.
  Previously a single missing phenotype value silently set every feature's
  coefficient to 0 and its p-value to 1.
* Spearman correlation again uses `cor.test()`, which supplies the exact
  small-sample p-value. The vectorised normal approximation that replaced it
  shifted p-values enough to change significance calls.
* `GetGroupMeans()` returns rows sorted by feature id again, so heatmap row
  order and exported tables match earlier releases. An unused factor level in
  the group column no longer produces an all-missing column.
* `GetFaithPD()`, `GetBetaDiversity()`, and `GetAlphaDiversity()` gained a
  `use_fastcluster` argument, default `FALSE`. Tree building previously used
  `fastcluster::hclust()` whenever that package happened to be installed;
  because it breaks ties between equal distances differently from
  `stats::hclust()`, the same data could give different diversity values on
  different machines.
* `GetBetaDiversity(method = "CSCS")` checks the similarity diagonal with
  `Matrix::diag()` instead of expanding the whole matrix to dense.
* `GenerateHeatmapInputs()` stops with an explanatory error above 10,000
  features rather than allocating a very large dense distance matrix.
* Tutorial 2 now uses `AddChemSim()`; it previously built a `.dissim` slot that
  the similarity-based diversity functions could not read. Tutorial 3 calls
  `FeatureDendrogram_derep()` and its note describes which functions read which
  slot.
* `AddChemSim()`, `GetSimMat()`, `AddCustomSim()`, and the six `_derep`
  functions are exported. The documentation block for `GetFaithPD()` had lost
  its title and was being attached to an internal helper.

# eCOMET 0.0.0.9000
* Initial packaged version from custom functions
* Roxygen-formated

# eCOMET 0.0.0.9001
* Vignettes added
* NMDS, PCoA added

# eCOMET 0.0.0.9002
* Importing functions now deal with mzML, mzXML, and thermo raw files
* All plot generating functions now returns a list of the plot and raw data. Saving the outputs are optional.

