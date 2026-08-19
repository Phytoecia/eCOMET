# What each test actually proves

Every expectation in this suite carries one of five trust grades. The grade says
what a failure means, which is not the same for all of them. Read the grade
before deciding whether a red test is a bug in the code or a change that needs a
new snapshot.

| Grade | What it proves | What a failure means |
|---|---|---|
| `hand` | The value was derived on paper from a 3x3 toy dataset and hard-coded. | The arithmetic no longer matches the published definition. Almost certainly a bug in eCOMET. **Do not edit the expected value.** |
| `external` | eCOMET agrees with an independent implementation it does **not** call. | eCOMET and the reference disagree. Investigate eCOMET first. |
| `regression` | A bug recorded in `NEWS.md` is still fixed. The correct answer is known because the fix already shipped. | The fix was undone. **Do not edit the expected value.** |
| `contract` | A structural invariant holds: dimensions, names, slot synchronisation, round trips, absence of file side effects. | The object contract changed. Decide whether that was intended. |
| `frozen` | Nothing about correctness. The value is a recording of current behaviour. | The output **changed**, not that it is now wrong. Investigate the change; if the new value is right, update the snapshot deliberately. |

Why `frozen` is separated out: freezing today's output turns today's bugs into
the specification. The Rao Q recycling bug in `NEWS.md` is the concrete example.
Had it been snapshotted at submission, the later fix would have shown up as a
test failure and could have been reverted to make the suite green. Grading the
snapshots keeps that from happening silently.

## Files

| File | Grade | Covers |
|---|---|---|
| `test-diversity-hand.R` | `hand` | `GetHillNumbers` (q = 0/1/2), `GetFunctionalHillNumber`, `GetBetaDiversity(method = "CSCS")`, `CanopusLevelEnrichmentAnal` |
| `test-external-parity.R` | `external` | `GetHillNumbers` vs `vegan`, `PairwiseComp` vs `stats::t.test(var.equal = FALSE)` + `p.adjust` |
| `test-regression-news.R` | `regression` | the eight fixes in `NEWS.md` |
| `test-mmo-contract.R` | `contract` | `filter_mmo`, `ReplaceZero`, the four normalisations, `GetNormFeature`, `SaveMMO`/`LoadMMO`, `FeatureToID`/`IDToFeature`, `print.mmo`, `save_output = FALSE` file side effects |
| `test-integration.R` | `contract` | MZmine CSV to PCoA in one pass, plus the column names the importers depend on |
| `test-derep-smoke.R` | `smoke` | the six superseded `_derep` functions: structure only, never values |
| `test-api-smoke.R` | `smoke` | exported functions not reached elsewhere; asserts that every name in `NAMESPACE` is called somewhere in the suite |
| `test-golden.R` | `frozen` | everything else |

## The hand-checked values

`helper-mmo.R` defines a three-feature, three-sample toy object whose numbers
were chosen so that every quantity has a closed form:

```
proportions   p1 = (1/3, 1/3, 1/3)   p2 = (1/2, 1/4, 1/4)   p3 = (3/5, 1/5, 1/5)
similarity    S(f1,f2) = 1/2   S(f1,f3) = 0   S(f2,f3) = 1/4   diagonal 1
Rao Q         Q1 = 1/2   Q2 = 0.53125   Q3 = 0.42
```

From those: Hill numbers are `3`, `2^(3/2)`, `25/11` and so on; the functional
Hill number of the uniform sample is `Q / p^2 = 4.5` at every order q (all Hill
numbers coincide for a uniform distribution); CSCS beta diversity is
`1 - 0.5/0.53125`, `1 - 0.5/0.58`, `1 - 0.55/0.58`; and the CANOPUS enrichment
p-value is the hypergeometric tail `280/1001`. Tolerance is `1e-8`.

## What is deliberately not tested

**Wrapper mathematics.** A large part of eCOMET calls out to `vegan`,
`picante`, `GUniFrac`, `fgsea`, `ape`, `caret`, `pairwiseAdonis` or
`stats::fisher.test`. Re-checking those results against the same package they
came from proves nothing. What *is* eCOMET's own code is the glue around the
call: transposition, conversion to relative proportions, row and column
alignment, name preservation, NA handling. The Rao Q bug lived precisely in that
glue. Tests on wrapper functions are written as transformation and alignment
checks, not as mathematical verification.

| Function | Delegates to | Source line |
|---|---|---|
| `GetBetaDiversity(method = "bray" / "jaccard")` | `vegan::vegdist` | 5955, 5966 |
| `GetBetaDiversity(method = "Gen.Uni")` | `GUniFrac::GUniFrac` | 5940 |
| `permanova_stat` | `vegan::adonis2`, `pairwiseAdonis::pairwise.adonis2` | 2515-2517 |
| `GetFaithPD` | `picante::pd` | 4835 |
| `MSEA` | `fgsea::fgsea` | 3642 |
| `CanopusLevelEnrichmentAnal` | `stats::fisher.test` | 3270 |
| `PLSDAplot` | `caret::plsda` | 2920 |
| `NMDSplot` | `vegan::metaMDS` | 6027 |
| `PCoAplot` | `ape::pcoa` | 6085 |
| correlation family | `stats::cor.test` | 3816, 3833, 3932, 4088 |

eCOMET's own implementations, where hand and external checks do carry weight:
`GetHillNumbers`, `GetFunctionalHillNumber`, the CSCS branch of
`GetBetaDiversity`, `GetRichness`, `GetSpecializationIndex`, `RarefactionAUC`,
`sparse_single_phylo`, and the vectorised Welch t-test in `PairwiseComp`.

**Plot images.** Plot functions return a list of the plot and its data. Only the
data component is snapshotted. Pixel comparison (vdiffr and similar) is
platform dependent and would keep CI red for reasons unrelated to the package.

**`_derep` versus `.sim` agreement.** The six `_derep` functions are the older
dissimilarity-based API and are now on an explicit `lifecycle::deprecate_soft()`
path. Asserting agreement would pressure the corrected new API back towards the
old behaviour, so the comparison lives in
`inst/dev/sim_vs_derep_divergence.md` as documentation instead. Measured there:
every pair agrees to floating point **except** the functional Hill family, where
`GetFunctionalHillNumber_derep()` divides by Rao Q twice and returns wrong
values (`GetFunctionalHillNumber()` reproduces the chemodiv 0.3.1 reference to
within 5e-15). That is a defect in `_derep`, not a difference of convention.

## Determinism

Anything stochastic is pinned before it is frozen.

| Function | Source of variation | Handling |
|---|---|---|
| `permanova_stat` | `permutations = 5000` | `set.seed(20260819)`, `permutations = 199` |
| `RarefactionAUC` / rarefied output | bootstrap, permutation | explicit `seed`, reduced `n_boot` / `n_perm` |
| `NMDSplot` | `metaMDS(try = 50, trymax = 100)` | coordinates not frozen; shape and column names checked instead |
| `PLSDAplot` | `caret::plsda` | `set.seed(20260819)`, tolerance `1e-4` |
| `MSEA` | `fgsea` permutations | `set.seed(20260819)`; only the deterministic enrichment score is frozen, p-values are checked for being probabilities |
| `hclust` | tie breaking differs between `stats` and `fastcluster` | `use_fastcluster = FALSE` throughout; the equivalence test builds its own tie-free matrix, because DreaMS scores are rounded to three decimals and the fixture necessarily contains ties |

Snapshot tolerance is `1e-6`, not for cosmetic reasons: `AddChemSim()` stores
dense similarity matrices with `mode(x) <- "single"`, so everything downstream of
a similarity matrix carries roughly float32 precision.

## Known defects surfaced as skips

Eight tests are `skip()`ped rather than asserted at current behaviour. Freezing a
defect would make it the specification. Each skip names the function, the source
line and what the function should do.

| Function | Defect | Source |
|---|---|---|
| `filter_mmo()` | `feature_presence` and `imputed_feature_data` are routed through the row-only branch, so their sample columns are never subset. `GetBetaDiversity(method = "jaccard", filter_group = TRUE)` returns a matrix over all samples. | `:2062-2069`, fix belongs at `:2041` |
| `ReplaceZero()` | `method = c("one", "half_min")` is never resolved by `match.arg()`, so the documented default call fails with `the condition has length > 1`. | `:822` |
| `VolcanoPlot()` | `print(plot)` runs unconditionally and leaves `Rplots.pdf` in the working directory even with `save_output = FALSE`. | `:2765` |
| `GetPerformanceFeatureLMM()`, `ScreenFeaturePhenotypeCorrelation(model = "lmm")` | Both read `summary(fit)$coefficients[2, 5]` from an `lme4::lmer` fit, whose coefficient matrix has three columns. The fifth column exists only under `lmerTest`, which is not a dependency, so both fail on every input. | `:3843`, `:3997` |
| `FeaturePhenotypeCorrelation(model = "lmm")` | The formula is built from `$`-extracted columns with no `data =` argument, so lme4 cannot resolve the grouping factor. | `:3703` |
| `GetFunctionalHillNumber_derep()` | Computes `Pq <- ratio^q` with `ratio = p / Q`, dividing by Rao Q twice, so it returns `FD(Q) * Q^(-q/(1-q))`. `GetAlphaDiversity_derep(mode = "weighted")` inherits the wrong values. | `:4498` |
| `GetPerformanceFeature*()` | The fold-change column is joined with `match(results$feature, mmo$pairwise$feature)`; the left side holds feature IDs and the right side holds `mz_rt` keys, so the column comes back entirely `NA`. The key should be `mmo$pairwise$id`. | `:3964`, `:4026`, `:4132` |

## Fixtures

`tests/testthat/fixtures/` is built by `data-raw/make_test_fixture.R` from the
tutorial data in `inst/extdata/tutorials/`. Regenerate with

```r
Rscript data-raw/make_test_fixture.R
```

after any change to the importers; never edit the fixtures by hand. The slice is
30 annotated features (plus three unannotated ones) across 9 samples in 3
groups, about 145 KB in total against 82 MB of tutorial data. It preserves the
properties the tests rely on: every fixture feature carries a CANOPUS
annotation, seven NPC pathways and 21 NPC classes are represented, six features
contain zeros, and the metadata carries both naturally missing phenotypes
(`sl`, `le_rgr`) and a constructed one-value-missing phenotype (`mass_missing`).

`mini_feature_table.csv` keeps MZmine's column headers verbatim. That is what
makes `test-integration.R` a canary for an MZmine format change, so never
rewrite those names.
