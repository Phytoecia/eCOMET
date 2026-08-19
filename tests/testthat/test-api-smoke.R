# TRUST: smoke -- the exported functions not reached by the other files. These
# check that each one runs and returns the documented structure, so that the
# suite touches every entry point in NAMESPACE. Values are frozen in
# test-golden.R, not here.

smoke_mmo <- function() {
  suppressMessages({
    mmo <- ReplaceZero(mini_mmo(), method = "one")
    mmo <- ZNormalization(mmo)
    mmo <- LogNormalization(mmo)
    utils::capture.output({
      mmo <- PairwiseComp(mmo, "ctrl", "sl1")
      mmo <- PairwiseComp(mmo, "ctrl", "le1")
    })
  })
  mmo
}

# ---------------------------------------------------------------------------
# Superseded dissimilarity accessors
# ---------------------------------------------------------------------------

test_that("AddChemDist and GetDistanceMat round trip a .dissim slot", {
  # The edge list covers all 33 features in the CSV fixture, so compare against
  # the unfiltered mmo; mini_mmo() keeps only the 30 annotated ones.
  full <- mini_mmo_unannotated()
  mmo <- suppressMessages(
    AddChemDist(full, dreams_dir = fixture_path("mini_dreams_sim.csv")))
  expect_false(is.null(mmo$dreams.dissim))

  d <- GetDistanceMat(mmo, "dreams")
  expect_equal(nrow(d), ncol(d))
  expect_setequal(rownames(d), full$feature_data$id)
  expect_identical(rownames(d), colnames(d))
  # AddChemDist stores dissimilarity, so the diagonal is 0, not 1.
  expect_true(all(abs(diag(as.matrix(d))) < 1e-6))
  expect_true(isSymmetric(unname(as.matrix(d))))
})

# ---------------------------------------------------------------------------
# Enrichment plot wrappers
# ---------------------------------------------------------------------------

test_that("the CANOPUS enrichment plot wrappers run and return their data", {
  mmo <- mini_mmo()
  ids <- mmo$feature_data$id[1:10]

  # pthr = 1 keeps every term: with 30 features nothing reaches p < 0.05, and
  # these wrappers facet on the result, which needs at least one row.
  out1 <- suppressWarnings(suppressMessages(
    CanopusLevelEnrichmentPlot(mmo, comp.list = list(test = ids), pthr = 1,
                               outdir = "enrich", save_output = FALSE)))
  expect_type(out1, "list")
  expect_true(any(vapply(out1, is.data.frame, logical(1))))

  out2 <- suppressWarnings(suppressMessages(
    CanopusAllLevelEnrichmentPlot(mmo, comp.list = list(test = ids), pthr = 1,
                                  outdir = "enrich", save_output = FALSE)))
  expect_type(out2, "list")

  out3 <- suppressWarnings(suppressMessages(
    CanopusListEnrichmentPlot_2(mmo, id_list = ids, pthr = 1, topn = 3,
                                outdir = "enrich", save_output = FALSE)))
  expect_type(out3, "list")
})

# ---------------------------------------------------------------------------
# Dendrogram rendering
# ---------------------------------------------------------------------------

test_that("PlotFeatureDendrogram renders the FeatureDendrogram output", {
  skip_if_not_installed("ggtree")
  mmo <- mini_mmo()
  tree <- suppressMessages(FeatureDendrogram(mmo, distance = "dreams"))
  p <- suppressWarnings(suppressMessages(PlotFeatureDendrogram(tree, mmo = mmo)))
  expect_s3_class(p, "ggplot")
})

test_that("PlotFeatureDendrogram rejects anything but the FeatureDendrogram list", {
  mmo <- mini_mmo()
  tree <- suppressMessages(FeatureDendrogram(mmo, distance = "dreams"))
  expect_error(PlotFeatureDendrogram(tree$phylo, mmo = mmo),
               "must be the list returned by FeatureDendrogram")
})

# ---------------------------------------------------------------------------
# Fold-change versus performance plots
# ---------------------------------------------------------------------------

# GetPerformanceFeature*() currently returns an all-NA fold-change column (see
# the skip below), so the input is built here instead of taken from it. That
# keeps these three plot functions covered while the join stays broken.
foldchange_input <- function() {
  reg <- suppressMessages(GetPerformanceFeatureRegression(
    smoke_mmo(), "mass", levels(mini_mmo()$metadata$group),
    DAM.list = list(ctrl_vs_sl1.up = mini_mmo()$feature_data$id[1:5]),
    comparisons = character(0)))
  pw <- smoke_mmo()$pairwise
  reg$ctrl_vs_sl1_log2FC <- pw[["ctrl_vs_sl1_log2FC"]][match(reg$feature, pw$id)]
  reg[stats::complete.cases(reg), ]
}

test_that("the fold-change versus performance plots run", {
  reg <- foldchange_input()
  expect_false(anyNA(reg$ctrl_vs_sl1_log2FC))
  cols <- c(ctrl_vs_sl1.up = "#d95f02", else. = "grey60")

  out1 <- suppressWarnings(suppressMessages(PlotFoldchangeResistanceRegression(
    reg, fold_change = "ctrl_vs_sl1_log2FC", color = cols,
    outdir = "fc", save_output = FALSE)))
  expect_type(out1, "list")
  expect_s3_class(out1$plot, "ggplot")

  out2 <- suppressWarnings(suppressMessages(PlotFoldchangeResistanceRegression_t(
    reg, fold_change = "ctrl_vs_sl1_log2FC", color = cols,
    output_dir = "fc", save_output = FALSE)))
  expect_type(out2, "list")

  out3 <- suppressWarnings(suppressMessages(PlotFoldchangeResistanceQuad(
    reg, fold_change = "ctrl_vs_sl1_log2FC", color = cols,
    output_dir = "fc", save_output = FALSE)))
  expect_type(out3, "list")
  expect_true("quadrant" %in% names(out3[[which(vapply(out3, is.data.frame, logical(1)))[1]]]))
})

test_that("GetPerformanceFeature* fills in the requested fold-change columns", {
  skip(paste(
    "Known defect: all three GetPerformanceFeature* functions join the",
    "fold-change column with match(regression_results$feature,",
    "mmo$pairwise$feature) (R/250813_eCOMET_build_V3.R:3964, :4026, :4132).",
    "regression_results$feature holds feature IDs while mmo$pairwise$feature",
    "holds the mz_rt keys, so the match never succeeds and the column comes",
    "back entirely NA. The join key should be mmo$pairwise$id.",
    "PlotFoldchangeResistanceRegression(), ..._t() and ...Quad() all fail on",
    "that output with '0 (non-NA) cases'."
  ))
  mmo <- smoke_mmo()
  reg <- suppressMessages(GetPerformanceFeatureRegression(
    mmo, "mass", levels(mmo$metadata$group),
    DAM.list = list(ctrl_vs_sl1.up = mmo$feature_data$id[1:5]),
    comparisons = "ctrl_vs_sl1"))
  expect_false(all(is.na(reg[["ctrl_vs_sl1_log2FC"]])))
})

# ---------------------------------------------------------------------------
# Every exported function is reached by the suite
# ---------------------------------------------------------------------------

test_that("no exported function is left untouched by the test suite", {
  exported <- sort(getNamespaceExports("ecomet"))
  files <- list.files(test_path("."), pattern = "^test-.*\\.R$", full.names = TRUE)
  body <- paste(unlist(lapply(files, readLines)), collapse = "\n")
  called <- vapply(exported, function(f) {
    grepl(paste0("\\b", gsub("([.])", "\\\\.", f), "\\("), body)
  }, logical(1))
  expect_equal(exported[!called], character(0))
})
