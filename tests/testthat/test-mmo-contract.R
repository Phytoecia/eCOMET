# TRUST: contract -- these tests assert structural invariants of the mmo object
# (dimensions, names, slot synchronisation, round trips, file side effects).
# They record no numeric results, so they stay valid across any change that
# preserves the object contract.
#
# Three known defects are surfaced as skip()s rather than being asserted at
# their current (wrong) behaviour. Freezing them would turn a bug into the
# specification. Each skip names the function and what it should do.

prepared_mmo <- function() {
  suppressMessages({
    mmo <- ReplaceZero(mini_mmo(), method = "one")
    mmo <- ZNormalization(mmo)
    mmo <- LogNormalization(mmo)
    mmo <- MeancenterNormalization(mmo)
  })
  mmo
}

# ---------------------------------------------------------------------------
# filter_mmo
# ---------------------------------------------------------------------------

test_that("filter_mmo(id_list) keeps every feature-indexed slot in step", {
  mmo <- mini_mmo()
  keep <- mmo$feature_data$id[1:10]
  sub <- suppressMessages(filter_mmo(mmo, id_list = keep, drop_empty_feat = FALSE))

  expect_setequal(sub$feature_data$id, keep)
  expect_setequal(sub$feature_info$id, keep)
  expect_setequal(sub$sirius_annot$id, keep)
  expect_setequal(sub$pairwise$id, keep)
  expect_setequal(sub$feature_presence$id, keep)

  sim <- GetSimMat(sub, "dreams")
  expect_equal(dim(sim), c(10L, 10L))
  expect_setequal(rownames(sim), keep)
  expect_identical(rownames(sim), colnames(sim))
})

test_that("filter_mmo(group_list) subsets metadata and feature_data columns", {
  mmo <- mini_mmo()
  sub <- suppressMessages(filter_mmo(mmo, group_list = "ctrl"))

  expect_equal(nrow(sub$metadata), 3L)
  expect_setequal(as.character(sub$metadata$group), "ctrl")
  expect_setequal(setdiff(names(sub$feature_data), c("id", "feature")),
                  sub$metadata$sample)
})

test_that("filter_mmo(group_list) also subsets feature_presence and imputed data", {
  skip(paste(
    "Known defect: filter_mmo() routes feature_presence and",
    "imputed_feature_data through the row-only branch",
    "(R/250813_eCOMET_build_V3.R:2062-2069) instead of",
    "filter_feature_sample_table(), so their sample columns are never subset.",
    "GetBetaDiversity(method='jaccard', filter_group=TRUE) therefore returns a",
    "matrix over ALL samples. Re-enable once the slots are added to the",
    "filter_feature_sample_table() loop at :2041."
  ))
  mmo <- suppressMessages(ReplaceZero(mini_mmo(), method = "one"))
  sub <- suppressMessages(filter_mmo(mmo, group_list = "ctrl"))
  expect_equal(names(sub$feature_presence), names(sub$feature_data))
  expect_equal(names(sub$imputed_feature_data), names(sub$feature_data))
})

test_that("GetBetaDiversity(jaccard) honours filter_group", {
  skip("Known defect: see the filter_mmo() skip above. Currently returns 9x9 for a 6-sample request.")
  beta <- suppressWarnings(suppressMessages(
    GetBetaDiversity(mini_mmo(), method = "jaccard", normalization = "PA",
                     filter_group = TRUE, group_list = c("ctrl", "sl1"))
  ))
  expect_equal(dim(beta), c(6L, 6L))
})

# ---------------------------------------------------------------------------
# ReplaceZero
# ---------------------------------------------------------------------------

test_that("ReplaceZero removes every zero and preserves the table shape", {
  mmo <- mini_mmo()
  for (method in c("one", "half_min")) {
    out <- suppressMessages(ReplaceZero(mmo, method = method))$imputed_feature_data
    mat <- as.matrix(out[, -(1:2), drop = FALSE])
    expect_false(any(mat == 0), info = method)
    expect_false(anyNA(mat), info = method)
    expect_equal(dim(out), dim(mmo$feature_data), info = method)
    expect_identical(names(out), names(mmo$feature_data), info = method)
    expect_identical(out$id, mmo$feature_data$id, info = method)
  }
})

test_that("ReplaceZero(method = 'one') writes exactly 1 where the input was 0", {
  mmo <- mini_mmo()
  before <- as.matrix(mmo$feature_data[, -(1:2), drop = FALSE])
  after <- as.matrix(
    suppressMessages(ReplaceZero(mmo, method = "one"))$imputed_feature_data[, -(1:2), drop = FALSE]
  )
  expect_true(any(before == 0))            # the fixture must exercise this path
  expect_true(all(after[before == 0] == 1))
  expect_equal(after[before != 0], before[before != 0])
})

test_that("ReplaceZero works without naming the method", {
  skip(paste(
    "Known defect: ReplaceZero(mmo) fails with 'the condition has length > 1'.",
    "method = c('one', 'half_min') is never resolved by match.arg()",
    "(R/250813_eCOMET_build_V3.R:822), so `if (method == 'one')` receives a",
    "length-2 condition, an error since R 4.2. Every tutorial passes method",
    "explicitly, which is why it has gone unnoticed."
  ))
  expect_no_error(suppressMessages(ReplaceZero(mini_mmo())))
})

# ---------------------------------------------------------------------------
# Normalisations
# ---------------------------------------------------------------------------

test_that("normalisations preserve dimensions and column names", {
  mmo <- mini_mmo()
  ref <- mmo$feature_data
  cases <- list(
    log          = function(m) suppressMessages(LogNormalization(m))$log,
    meancentered = function(m) suppressMessages(MeancenterNormalization(m))$meancentered,
    zscore       = function(m) suppressMessages(ZNormalization(m))$zscore
  )
  for (nm in names(cases)) {
    out <- cases[[nm]](mmo)
    expect_equal(dim(out), dim(ref), info = nm)
    expect_identical(names(out), names(ref), info = nm)
    expect_identical(out$id, ref$id, info = nm)
  }
})

test_that("MassNormalization rewrites feature_data in place with the same shape", {
  mmo <- mini_mmo()
  out <- suppressMessages(MassNormalization(mmo))$feature_data
  expect_equal(dim(out), dim(mmo$feature_data))
  expect_identical(names(out), names(mmo$feature_data))
  expect_identical(out$id, mmo$feature_data$id)
})

test_that("GetNormFeature returns the slot each normalisation wrote", {
  mmo <- prepared_mmo()
  expect_identical(GetNormFeature(mmo, "None"), mmo$feature_data)
  expect_identical(GetNormFeature(mmo, "PA"), mmo$feature_presence)
  expect_identical(GetNormFeature(mmo, "Log"), mmo$log)
  expect_identical(GetNormFeature(mmo, "Z"), mmo$zscore)
  expect_identical(GetNormFeature(mmo, "Meancentered"), mmo$meancentered)
  expect_identical(GetNormFeature(mmo, "Imputed"), mmo$imputed_feature_data)
})

# ---------------------------------------------------------------------------
# Identifier round trips, printing, serialisation
# ---------------------------------------------------------------------------

test_that("FeatureToID and IDToFeature invert each other", {
  mmo <- mini_mmo()
  ids <- mmo$feature_data$id[1:5]
  expect_equal(as.character(FeatureToID(mmo, IDToFeature(mmo, ids))), ids)

  features <- mmo$feature_data$feature[1:5]
  expect_equal(as.character(IDToFeature(mmo, FeatureToID(mmo, features))), features)
})

test_that("SaveMMO and LoadMMO round trip the object", {
  mmo <- mini_mmo()
  withr::with_tempdir({
    SaveMMO(mmo, file = "mmo.rds", include_session = FALSE)
    back <- suppressMessages(LoadMMO("mmo.rds", verbose = FALSE))
    expect_equal(back, mmo)
  })
})

test_that("SaveMMO(include_session = TRUE) adds only the session attribute", {
  mmo <- mini_mmo()
  withr::with_tempdir({
    SaveMMO(mmo, file = "mmo.rds", include_session = TRUE)
    back <- suppressMessages(LoadMMO("mmo.rds", verbose = FALSE))
    expect_true(!is.null(attr(back, "saved_session_info")))
    attr(back, "saved_session_info") <- NULL
    expect_equal(back, mmo)
  })
})

test_that("print.mmo runs and mentions the object size", {
  expect_output(print(mini_mmo()), "feature", ignore.case = TRUE)
})

# ---------------------------------------------------------------------------
# save_output = FALSE must not touch the filesystem
# ---------------------------------------------------------------------------

no_files_written <- function(code) {
  withr::with_tempdir({
    before <- list.files(".", recursive = TRUE, all.files = TRUE, no.. = TRUE)
    # invisible() matters: capture.output() prints the value of its expression,
    # and printing a returned ggplot would render it and open a device here.
    suppressWarnings(suppressMessages(utils::capture.output(invisible(force(code)))))
    after <- list.files(".", recursive = TRUE, all.files = TRUE, no.. = TRUE)
    setdiff(after, before)
  })
}

test_that("plot functions write nothing when save_output = FALSE", {
  mmo <- prepared_mmo()
  beta <- suppressMessages(GetBetaDiversity(mmo, method = "bray", normalization = "None"))
  ids <- mmo$feature_data$id

  expect_equal(no_files_written(PCAplot(mmo, outdir = "PCA", save_output = FALSE)),
               character(0))
  expect_equal(no_files_written(PCoAplot(mmo, betadiv = beta, outdir = "pcoa", save_output = FALSE)),
               character(0))
  expect_equal(no_files_written({
    set.seed(1)
    NMDSplot(mmo, betadiv = beta, outdir = "nmds", save_output = FALSE)
  }), character(0))
  expect_equal(no_files_written({
    set.seed(1)
    PLSDAplot(mmo, color = c(ctrl = "blue", sl1 = "red", le1 = "green"),
              outdir = "plsda", save_output = FALSE)
  }), character(0))
  expect_equal(no_files_written(
    PlotNPCStackedBar(mmo, group_col = "group", outdir = "npc", save_output = FALSE)),
    character(0))
  expect_equal(no_files_written(
    AnovaBarPlot(mmo, ID_list = ids[1:2], outdir = "anova", save_output = FALSE)),
    character(0))
  expect_equal(no_files_written(
    CanopusListEnrichmentPlot(mmo, id_list = ids[1:8], pthr = 1,
                              outdir = "enrich", save_output = FALSE)),
    character(0))
  expect_equal(no_files_written({
    set.seed(1)
    MSEA(mmo, feature_id = ids, feature_score = seq_along(ids),
         term_level = "NPC_pathway", outdir = "msea", save_output = FALSE)
  }), character(0))
})

test_that("VolcanoPlot writes nothing when save_output = FALSE", {
  skip(paste(
    "Known defect: VolcanoPlot() calls print(plot) unconditionally",
    "(R/250813_eCOMET_build_V3.R:2765), which opens the default graphics",
    "device and leaves Rplots.pdf in the working directory in a",
    "non-interactive session even when save_output = FALSE."
  ))
  mmo <- prepared_mmo()
  suppressMessages(utils::capture.output(mmo <- PairwiseComp(mmo, "ctrl", "sl1")))
  expect_equal(
    no_files_written(VolcanoPlot(mmo, comp = "ctrl_vs_sl1", outdir = "v.png",
                                 save_output = FALSE)),
    character(0)
  )
})
