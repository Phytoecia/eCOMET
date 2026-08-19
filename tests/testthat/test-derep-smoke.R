# TRUST: smoke -- the six _derep functions are the superseded dissimilarity-based
# API (_pkgdown.yml: "Dissimilarity-based versions (superseded)") and are now on
# a lifecycle::deprecate_soft() path.
#
# This file checks only that they RUN and return the documented structure:
# class, dimensions, dimnames, column names. It deliberately does NOT freeze
# their values, because doing so would pressure the corrected new API back
# towards the old behaviour. The measured value comparison lives in
# inst/dev/sim_vs_derep_divergence.md as documentation.
#
# One of the six is known to be wrong rather than merely old; see the skipped
# test at the bottom of this file.

derep_mmo <- function() {
  mmo <- mini_mmo()
  dis <- 1 - as.matrix(GetSimMat(mmo, "dreams"))
  mmo <- suppressMessages(AddCustomDist(mmo, dis, name = "chem"))
  suppressMessages(ReplaceZero(mmo, method = "one"))
}

derep_dist <- function(mmo = derep_mmo()) mmo$chem.dissim

test_that("GetFunctionalHillNumber_derep returns one row per sample", {
  withr::local_options(lifecycle_verbosity = "quiet")
  mmo <- derep_mmo()
  for (q in c(0, 1, 2)) {
    out <- GetFunctionalHillNumber_derep(mmo$feature_data, mmo$metadata,
                                         derep_dist(mmo), q = q)
    expect_s3_class(out, "data.frame")
    expect_named(out, c("sample", "group", "hill_number", "value"))
    expect_equal(nrow(out), nrow(mmo$metadata), info = paste("q =", q))
    expect_setequal(out$sample, mmo$metadata$sample)
    expect_type(out$hill_number, "double")
  }
})

test_that("GetFaithPD_derep returns picante's pd columns plus sample and group", {
  withr::local_options(lifecycle_verbosity = "quiet")
  mmo <- derep_mmo()
  out <- GetFaithPD_derep(mmo$feature_data, mmo$metadata, derep_dist(mmo))
  expect_s3_class(out, "data.frame")
  expect_true(all(c("PD", "SR", "sample", "group", "value") %in% names(out)))
  expect_equal(nrow(out), nrow(mmo$metadata))
  expect_setequal(out$sample, mmo$metadata$sample)
})

test_that("GetAlphaDiversity_derep runs for each mode and returns per-sample rows", {
  withr::local_options(lifecycle_verbosity = "quiet")
  mmo <- derep_mmo()
  for (mode in c("richness", "unweighted", "weighted")) {
    out <- suppressMessages(
      GetAlphaDiversity_derep(mmo, mode = mode, distance = "chem",
                              normalization = "None", q = 1)
    )
    expect_s3_class(out, "data.frame")
    expect_true(all(c("sample", "group", "value") %in% names(out)), info = mode)
    expect_equal(nrow(out), nrow(mmo$metadata), info = mode)
    expect_setequal(out$sample, mmo$metadata$sample)
  }
})

test_that("GetBetaDiversity_derep returns square sample matrices", {
  withr::local_options(lifecycle_verbosity = "quiet")
  mmo <- derep_mmo()
  samples <- mmo$metadata$sample

  bray <- suppressMessages(
    GetBetaDiversity_derep(mmo, method = "bray", normalization = "None"))
  expect_true(is.matrix(bray))
  expect_equal(dim(bray), c(length(samples), length(samples)))
  expect_setequal(rownames(bray), samples)
  expect_identical(rownames(bray), colnames(bray))

  cscs <- suppressMessages(
    GetBetaDiversity_derep(mmo, method = "CSCS", distance = "chem",
                           normalization = "None"))
  expect_true(is.matrix(cscs))
  expect_equal(dim(cscs), c(length(samples), length(samples)))
  expect_identical(rownames(cscs), colnames(cscs))

  guni <- suppressMessages(
    GetBetaDiversity_derep(mmo, method = "Gen.Uni", distance = "chem",
                           normalization = "None"))
  expect_type(guni, "list")
  expect_named(guni, c("d_0", "d_0.5", "d_1"))
  expect_equal(dim(guni$d_0.5), c(length(samples), length(samples)))
})

test_that("GenerateHeatmapInputs_derep returns the four heatmap components", {
  withr::local_options(lifecycle_verbosity = "quiet")
  mmo <- derep_mmo()
  out <- suppressMessages(
    GenerateHeatmapInputs_derep(mmo, distance = "chem", normalization = "None"))
  expect_named(out, c("FC_matrix", "dist_matrix", "row_label", "heatmap_data"))
  expect_true(is.matrix(out$FC_matrix))
  expect_equal(nrow(out$FC_matrix), nrow(mmo$feature_data))
  expect_s3_class(out$dist_matrix, "dist")
  expect_equal(attr(out$dist_matrix, "Size"), nrow(out$FC_matrix))
  expect_equal(out$row_label, rownames(out$FC_matrix))
})

test_that("FeatureDendrogram_derep returns a tree over the fixture features", {
  withr::local_options(lifecycle_verbosity = "quiet")
  mmo <- derep_mmo()
  out <- suppressMessages(
    FeatureDendrogram_derep(mmo, distance = "chem", save_newick = FALSE))
  tree <- if (inherits(out, "phylo")) out else out[[which(vapply(out, inherits, logical(1), "phylo"))[1]]]
  expect_s3_class(tree, "phylo")
  expect_setequal(tree$tip.label, mmo$feature_data$id)
})

test_that("the _derep functions are exported alongside the .sim entry points", {
  withr::local_options(lifecycle_verbosity = "quiet")
  derep <- c("GetAlphaDiversity_derep", "GetBetaDiversity_derep", "GetFaithPD_derep",
             "GetFunctionalHillNumber_derep", "FeatureDendrogram_derep",
             "GenerateHeatmapInputs_derep")
  exported <- getNamespaceExports("ecomet")
  expect_true(all(derep %in% exported))
  expect_true(all(sub("_derep$", "", derep) %in% exported))
})

test_that("each _derep function announces its own deprecation", {
  # The quiet option above hides the warning while the structure is checked.
  # Here it is the thing being checked: the six functions must be on an
  # explicit deprecation path, not silently superseded.
  mmo <- withr::with_options(list(lifecycle_verbosity = "quiet"), derep_mmo())
  dist <- mmo$chem.dissim

  expect_warning(
    GetFunctionalHillNumber_derep(mmo$feature_data, mmo$metadata, dist, q = 1),
    class = "lifecycle_warning_deprecated")
  expect_warning(
    GetFaithPD_derep(mmo$feature_data, mmo$metadata, dist),
    class = "lifecycle_warning_deprecated")
  expect_warning(
    suppressMessages(GetAlphaDiversity_derep(mmo, mode = "richness", distance = "chem")),
    class = "lifecycle_warning_deprecated")
  expect_warning(
    suppressMessages(GetBetaDiversity_derep(mmo, method = "bray")),
    class = "lifecycle_warning_deprecated")
  expect_warning(
    suppressMessages(GenerateHeatmapInputs_derep(mmo, distance = "chem")),
    class = "lifecycle_warning_deprecated")
  expect_warning(
    suppressMessages(FeatureDendrogram_derep(mmo, distance = "chem")),
    class = "lifecycle_warning_deprecated")
})

test_that("GetFunctionalHillNumber_derep agrees with the corrected implementation", {
  skip(paste(
    "Known defect: GetFunctionalHillNumber_derep() computes Pq <- ratio^q with",
    "ratio = p / Q (R/250813_eCOMET_build_V3.R:4498), dividing by Rao Q twice.",
    "It returns FD(Q) * Q^(-q/(1-q)), which is Q * FD(Q) at q = 1. Commit",
    "90298d2 corrected GetFunctionalHillNumber() and left this one alone;",
    "GetFunctionalHillNumber() was since verified against the chemodiv 0.3.1",
    "reference implementation to within 5e-15. GetAlphaDiversity_derep() still",
    "routes through this function, so GetAlphaDiversity_derep(mode = 'weighted')",
    "returns wrong numbers. Measured divergence:",
    "inst/dev/sim_vs_derep_divergence.md."
  ))
  mmo <- derep_mmo()
  for (q in c(0, 1, 2)) {
    new <- GetFunctionalHillNumber(mmo$feature_data, mmo$metadata,
                                   GetSimMat(mmo, "dreams"), q = q)$hill_number
    old <- GetFunctionalHillNumber_derep(mmo$feature_data, mmo$metadata,
                                         derep_dist(mmo), q = q)$hill_number
    expect_equal(old, new, tolerance = 1e-6, info = paste("q =", q))
  }
})
