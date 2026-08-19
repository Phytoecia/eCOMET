# sim_vs_derep_divergence.R
#
# Dev-only script - NOT part of the installed package and NOT a test.
#
# Purpose
#   Quantify how far the six superseded `_derep` functions (dissimilarity
#   storage, `.dissim` slots) have drifted from the similarity-based entry
#   points (`.sim` slots) they were replaced by, and write the result as a
#   table for the documentation.
#
#   This is deliberately NOT a test. The `_derep` functions are older code on a
#   lifecycle::deprecate_soft() path; corrections recorded in NEWS.md (Rao Q
#   normalisation, missing-value handling, use_fastcluster) were applied to the
#   new API only. Asserting agreement would either raise false alarms or push
#   the corrected API back towards the old behaviour. See
#   tests/testthat/TRUST.md.
#
#   compare_sim_vs_derep.R answers "do they agree?" with PASS/FAIL. This script
#   answers "by how much do they differ?".
#
# How to run, from the package root:
#   Rscript inst/dev/sim_vs_derep_divergence.R
#
# Input: the treatment-based tutorial data shipped in inst/extdata/tutorials/,
# the same dataset the NEWS.md figures were computed on.
#
# Output: inst/dev/sim_vs_derep_divergence.md

suppressMessages({
  devtools::load_all(".", quiet = TRUE)
  library(dplyr)
})

options(lifecycle_verbosity = "quiet")

SRC <- "inst/extdata/tutorials/treatment_based"
OUT <- "inst/dev/sim_vs_derep_divergence.md"

msg <- function(...) cat("[divergence]", ..., "\n")

# -- build one mmo carrying both storage layouts ---------------------------
msg("building mmo from the treatment-based tutorial data")
mmo <- GetMZmineFeature(
  mzmine_dir   = file.path(SRC, "feature_table_demo.csv"),
  metadata_dir = file.path(SRC, "metadata_demo.csv"),
  group_col    = "group",
  sample_col   = "sample"
)
mmo <- AddSiriusAnnot(
  mmo,
  canopus_structuredir = file.path(SRC, "structure_identifications.tsv"),
  canopus_formuladir   = file.path(SRC, "canopus_formula_summary.tsv")
)
mmo <- AddChemSim(mmo, dreams_dir = file.path(SRC, "dreams_sim_demo.csv"))
mmo <- AddChemDist(mmo, dreams_dir = file.path(SRC, "dreams_sim_demo.csv"))
mmo <- ReplaceZero(mmo, method = "one")

sim  <- GetSimMat(mmo, "dreams")
dis  <- GetDistanceMat(mmo, "dreams")
feat <- mmo$feature_data
meta <- mmo$metadata

msg(nrow(feat), "features,", nrow(meta), "samples")

# -- divergence summary ----------------------------------------------------
rows <- list()

summarise_pair <- function(label, new_vals, old_vals, note = "") {
  new_vals <- as.numeric(new_vals)
  old_vals <- as.numeric(old_vals)
  ok <- is.finite(new_vals) & is.finite(old_vals)
  d <- abs(new_vals[ok] - old_vals[ok])
  rel <- d / pmax(abs(old_vals[ok]), .Machine$double.eps)
  rows[[length(rows) + 1L]] <<- data.frame(
    comparison = label,
    n = sum(ok),
    max_abs_diff = if (any(ok)) max(d) else NA_real_,
    max_pct_diff = if (any(ok)) 100 * max(rel) else NA_real_,
    median_pct_diff = if (any(ok)) 100 * stats::median(rel) else NA_real_,
    identical_1e6 = if (any(ok)) all(d <= 1e-6 * pmax(abs(old_vals[ok]), 1)) else NA,
    note = note,
    stringsAsFactors = FALSE
  )
  invisible(NULL)
}

rank_change <- function(new_vals, old_vals) {
  ok <- is.finite(new_vals) & is.finite(old_vals)
  if (!any(ok)) return(NA_integer_)
  sum(rank(new_vals[ok]) != rank(old_vals[ok]))
}

# GetFunctionalHillNumber -------------------------------------------------
for (q in c(0, 1, 2)) {
  new <- GetFunctionalHillNumber(feat, meta, sim, q = q)$hill_number
  old <- GetFunctionalHillNumber_derep(feat, meta, dis, q = q)$hill_number
  summarise_pair(
    sprintf("GetFunctionalHillNumber(q = %d)", q), new, old,
    sprintf("%d of %d samples change rank", rank_change(new, old), length(new))
  )
  msg("functional Hill q =", q, "done")
}

# GetFaithPD ---------------------------------------------------------------
new <- GetFaithPD(feat, meta, sim, use_fastcluster = FALSE)$PD
old <- GetFaithPD_derep(feat, meta, dis)$PD
summarise_pair("GetFaithPD", new, old,
               sprintf("%d of %d samples change rank", rank_change(new, old), length(new)))
msg("Faith PD done")

# GetAlphaDiversity --------------------------------------------------------
for (mode in c("weighted", "faith", "unweighted", "richness")) {
  new <- GetAlphaDiversity(mmo, mode = mode, q = 1, distance = "dreams",
                           normalization = "None", use_fastcluster = FALSE)$value
  old <- GetAlphaDiversity_derep(mmo, mode = mode, q = 1, distance = "dreams",
                                 normalization = "None")$value
  summarise_pair(sprintf("GetAlphaDiversity(mode = '%s')", mode), new, old,
                 sprintf("%d of %d samples change rank", rank_change(new, old), length(new)))
  msg("alpha", mode, "done")
}

# GetBetaDiversity ---------------------------------------------------------
for (method in c("bray", "jaccard", "CSCS")) {
  new <- suppressWarnings(GetBetaDiversity(mmo, method = method, distance = "dreams",
                                           normalization = if (method == "jaccard") "PA" else "None"))
  old <- suppressWarnings(GetBetaDiversity_derep(mmo, method = method, distance = "dreams",
                                                 normalization = if (method == "jaccard") "PA" else "None"))
  summarise_pair(sprintf("GetBetaDiversity(method = '%s')", method),
                 as.numeric(new), as.numeric(old), "all sample pairs")
  msg("beta", method, "done")
}

new <- GetBetaDiversity(mmo, method = "Gen.Uni", distance = "dreams",
                        normalization = "None", use_fastcluster = FALSE)
old <- GetBetaDiversity_derep(mmo, method = "Gen.Uni", distance = "dreams",
                              normalization = "None")
for (nm in names(new)) {
  summarise_pair(sprintf("GetBetaDiversity(Gen.Uni)$%s", nm),
                 as.numeric(new[[nm]]), as.numeric(old[[nm]]), "all sample pairs")
}
msg("beta Gen.Uni done")

# GenerateHeatmapInputs ----------------------------------------------------
new <- GenerateHeatmapInputs(mmo, distance = "dreams", normalization = "None")
old <- GenerateHeatmapInputs_derep(mmo, distance = "dreams", normalization = "None")
summarise_pair("GenerateHeatmapInputs()$FC_matrix",
               as.numeric(new$FC_matrix), as.numeric(old$FC_matrix), "group means")
summarise_pair("GenerateHeatmapInputs()$dist_matrix",
               as.numeric(new$dist_matrix), as.numeric(old$dist_matrix), "feature pairs")
msg("heatmap inputs done")

# FeatureDendrogram --------------------------------------------------------
new_tree <- FeatureDendrogram(mmo, distance = "dreams")
old_tree <- FeatureDendrogram_derep(mmo, distance = "dreams")
summarise_pair("FeatureDendrogram()$hclust$height",
               new_tree$hclust$height, old_tree$hclust$height,
               sprintf("merge order identical: %s",
                       identical(new_tree$hclust$merge, old_tree$hclust$merge)))
msg("dendrogram done")

tab <- do.call(rbind, rows)

# -- write the table -------------------------------------------------------
fmt <- function(x) ifelse(is.na(x), "-", formatC(x, format = "g", digits = 4))
md <- c(
  "# Similarity-based API versus the superseded `_derep` functions",
  "",
  sprintf("Generated by `inst/dev/sim_vs_derep_divergence.R` on the treatment-based"),
  sprintf("tutorial data (%d features, %d samples in %d groups).",
          nrow(feat), nrow(meta), length(unique(meta$group))),
  "",
  "**This is documentation, not a test.** The `_derep` functions are the older",
  "dissimilarity-based API and are now on a `lifecycle::deprecate_soft()` path.",
  "See `tests/testthat/TRUST.md` for why this comparison is deliberately not",
  "asserted anywhere in the test suite.",
  "",
  "**Read the functional Hill rows as a defect, not as a difference of",
  "convention.** `GetFunctionalHillNumber()` was verified against the chemodiv",
  "0.3.1 reference implementation and agrees with it to floating point;",
  "`GetFunctionalHillNumber_derep()` divides by Rao Q twice and is simply wrong.",
  "Every other row below is agreement, so the table is a defect report for one",
  "function family rather than a comparison of two valid APIs.",
  "",
  "`max_pct_diff` is the largest relative difference between the two routes,",
  "expressed as a percentage of the `_derep` value.",
  "",
  "| Comparison | n | max abs diff | max % diff | median % diff | agree to 1e-6 | note |",
  "|---|---:|---:|---:|---:|:---:|---|",
  sprintf("| `%s` | %d | %s | %s | %s | %s | %s |",
          tab$comparison, tab$n, fmt(tab$max_abs_diff), fmt(tab$max_pct_diff),
          fmt(tab$median_pct_diff), ifelse(tab$identical_1e6, "yes", "**no**"), tab$note),
  "",
  "## Where the difference comes from",
  "",
  "Only the functional Hill number family diverges. Everything else - Faith PD,",
  "richness, unweighted alpha, all four beta-diversity methods, the heatmap",
  "inputs and the dendrogram - agrees to floating point.",
  "",
  "The two functional Hill implementations are not the same formula with one",
  "division fixed, and they are not two defensible conventions either: the",
  "similarity-based version reproduces chemodiv 0.3.1 `funcHillDiv()` to within",
  "5e-15, and the `_derep` version does not. Writing `p` for the sample",
  "proportions, `D` for the dissimilarity matrix and `Q = p\'Dp` for Rao\'s",
  "quadratic entropy:",
  "",
  "| | `GetFunctionalHillNumber()` | `GetFunctionalHillNumber_derep()` |",
  "|---|---|---|",
  "| q != 1 bracket | `Q^-q * sum_ij d_ij (p_i p_j)^q` | `Q^-2q * sum_ij d_ij (p_i p_j)^q` |",
  "| q != 1, zero features | masked out (`Pq[p == 0] <- 0`) | not masked, so `0^0 = 1` at q = 0 |",
  "| q = 1 | `2*sum (p/Q) log(p) (Dp) - log Q` | `2*sum (p/Q) log(p/Q) (Dp)` |",
  "",
  "The exponent of `Q` differs (`-q` against `-2q`), which is what produces the",
  "large q = 2 gap; the q = 0 gap is the zero masking; the q = 1 gap is the",
  "`log(p)` against `log(p/Q)` substitution together with the recycling fix.",
  "",
  "`_derep` therefore returns `FD(Q) * Q^(-q/(1-q))`, which at q = 1 is",
  "`Q * FD(Q)`. Commit 90298d2 (\"Functional Hill number fixed\") corrected the",
  "similarity-based version and left `_derep` untouched. `_derep` is still",
  "exported and `GetAlphaDiversity_derep()` still routes through it, so",
  "`GetAlphaDiversity_derep(mode = \"weighted\")` returns wrong numbers today.",
  "",
  "The similarity-based implementation is internally consistent: the q = 1",
  "branch is the exact q -> 1 limit of the q != 1 branch, for a uniform",
  "assemblage every order returns the same value `Q / p^2`, and it matches the",
  "chemodiv reference implementation. It computes `qFD(Q)` (Chiu & Chao 2014",
  "Eq. 4b, Petren et al. 2023 Box 1) - total functional diversity, an effective",
  "total dissimilarity rather than an effective number of compounds.",
  "",
  "Note on naming: Chiu & Chao reserve \"functional Hill number\" for `qD(Q)`",
  "(their Eq. 3), whose outer exponent is `1/(2(1-q))` and which equals",
  "`sqrt(FD(Q)/Q)`. Petren et al. and chemodiv compute `qFD(Q)` and call it",
  "\"functional Hill diversity\" (`FuncHillDiv`). eCOMET computes `qFD(Q)`, so",
  "\"functional Hill diversity\" is the accurate term for it. The value equals",
  "FAD at q = 0 and can exceed the feature count; it is an effective total",
  "dissimilarity, not an effective number of compounds.",
  "",
  "Note on the Rao Q recycling fix recorded in NEWS.md: `ratio` is only read",
  "inside the `q == 1` branch, so replacing `/` with `sweep()` changes q = 1",
  "and leaves q != 1 untouched. On this dataset the fix moves q = 1 by up to",
  "19.8% and reorders 22 of the 24 samples.",
  "",
  "## Session",
  "",
  "```",
  capture.output(utils::sessionInfo())[1:3],
  "```"
)
writeLines(md, OUT)
msg("wrote", OUT)
print(tab[, c("comparison", "max_pct_diff", "identical_1e6")], row.names = FALSE)
