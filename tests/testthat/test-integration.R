# TRUST: contract -- one pass through the whole import-to-ordination pipeline,
# starting from a CSV that carries MZmine's own column headers verbatim.
#
# This is the file that answers "the pipeline is likely to break with each new
# version of mzmine or SIRIUS": if either tool renames or drops a column the
# importers read, this test fails and names the step. Each stage asserts what it
# produced, so a failure localises the break instead of only reporting that the
# run died.
#
# fixtures/mini_feature_table.csv keeps the MZmine header strings exactly as
# exported (datafile:<sample>.mzML:area and friends). Never rewrite them: doing
# so would make the canary stop singing.

pipeline_paths <- function() {
  list(
    features = fixture_path("mini_feature_table.csv"),
    metadata = fixture_path("mini_metadata.csv"),
    canopus  = fixture_path("mini_canopus.tsv"),
    structure = fixture_path("mini_structure.tsv"),
    dreams   = fixture_path("mini_dreams_sim.csv")
  )
}

test_that("the fixture still carries the MZmine columns the importer reads", {
  hdr <- names(data.table::fread(pipeline_paths()$features, nrows = 0,
                                 check.names = FALSE))
  # Feature-level columns GetMZmineFeature() requires or requests by default.
  expect_true(all(c("id", "mz", "rt") %in% hdr))
  expect_true(all(c("rt_range:min", "rt_range:max",
                    "mz_range:min", "mz_range:max",
                    "feature_group",
                    "ion_identities:iin_id",
                    "ion_identities:ion_identities") %in% hdr))
  # Per-sample area columns, the ones fuzzy-matched against metadata.
  area <- grep("^datafile:.*:area$", hdr, value = TRUE)
  expect_equal(length(area), 9L)
})

test_that("the SIRIUS fixtures still carry the columns AddSiriusAnnot joins on", {
  p <- pipeline_paths()
  can <- names(data.table::fread(p$canopus, nrows = 0, check.names = FALSE))
  str_id <- names(data.table::fread(p$structure, nrows = 0, check.names = FALSE))
  expect_true("mappingFeatureId" %in% can)
  expect_true("mappingFeatureId" %in% str_id)
  expect_true(all(c("NPC#pathway", "NPC#superclass", "NPC#class") %in% can))
  expect_true(all(c("InChIkey2D", "smiles", "name") %in% str_id))
})

test_that("the DreaMS fixture is still the five-column edge list AddChemSim reads", {
  edges <- data.table::fread(pipeline_paths()$dreams, nrows = 5, check.names = FALSE)
  expect_equal(ncol(edges), 5L)
  expect_equal(names(edges)[1:2], c("ID1", "ID2"))
})

test_that("MZmine CSV to PCoA runs end to end", {
  p <- pipeline_paths()

  # 1. import ------------------------------------------------------------
  mmo <- suppressMessages(GetMZmineFeature(
    mzmine_dir   = p$features,
    metadata_dir = p$metadata,
    group_col    = "group",
    sample_col   = "sample"
  ))
  expect_s3_class(mmo, "mmo")
  expect_equal(nrow(mmo$feature_data), 33L)
  expect_equal(nrow(mmo$metadata), 9L)
  # The metadata file starts with a UTF-8 BOM; fread strips it, so the sample
  # column must still be addressable under its plain name.
  expect_true("sample" %in% names(mmo$metadata))
  expect_setequal(setdiff(names(mmo$feature_data), c("id", "feature")),
                  mmo$metadata$sample)
  expect_equal(nlevels(droplevels(mmo$metadata$group)), 3L)

  # 2. zero handling ------------------------------------------------------
  mmo <- suppressMessages(ReplaceZero(mmo, method = "one"))
  expect_false(any(as.matrix(mmo$imputed_feature_data[, -(1:2)]) == 0))

  # 3. SIRIUS annotation --------------------------------------------------
  mmo <- suppressMessages(AddSiriusAnnot(
    mmo,
    canopus_structuredir = p$structure,
    canopus_formuladir   = p$canopus
  ))
  expect_true("sirius_annot" %in% names(mmo))
  expect_equal(nrow(mmo$sirius_annot), nrow(mmo$feature_data))
  annotated <- sum(!is.na(mmo$sirius_annot[["NPC#pathway"]]))
  expect_equal(annotated, 30L)   # 30 annotated + 3 deliberately unannotated
  expect_gte(length(unique(na.omit(mmo$sirius_annot[["NPC#pathway"]]))), 3L)

  # 4. spectral similarity ------------------------------------------------
  mmo <- suppressMessages(AddChemSim(mmo, dreams_dir = p$dreams))
  sim <- GetSimMat(mmo, "dreams")
  expect_equal(nrow(sim), ncol(sim))
  expect_setequal(rownames(sim), mmo$feature_data$id)
  expect_true(all(abs(Matrix::diag(sim) - 1) < 1e-6))

  # 5. alpha diversity ----------------------------------------------------
  alpha <- suppressMessages(GetAlphaDiversity(
    mmo, mode = "weighted", q = 1, distance = "dreams", normalization = "None"))
  expect_s3_class(alpha, "data.frame")
  expect_equal(nrow(alpha), 9L)
  expect_setequal(alpha$sample, mmo$metadata$sample)
  expect_false(anyNA(alpha$value))

  # 6. beta diversity -----------------------------------------------------
  beta <- suppressMessages(GetBetaDiversity(
    mmo, method = "CSCS", distance = "dreams", normalization = "None"))
  expect_equal(dim(beta), c(9L, 9L))
  expect_true(isSymmetric(beta))
  expect_equal(unname(diag(beta)), rep(0, 9), tolerance = 1e-8)

  # 7. ordination ---------------------------------------------------------
  pcoa <- suppressWarnings(suppressMessages(
    PCoAplot(mmo, betadiv = beta, outdir = "pcoa", save_output = FALSE)))
  expect_named(pcoa, c("plot", "df", "permanova"))
  expect_s3_class(pcoa$plot, "ggplot")
  expect_equal(nrow(pcoa$df), 9L)
  expect_true(all(c("PCoA1", "PCoA2", "group") %in% names(pcoa$df)))
})

test_that("GetMZmineFeature reports a metadata sample with no matching area column", {
  p <- pipeline_paths()
  withr::with_tempdir({
    meta <- data.table::fread(p$metadata, check.names = FALSE, data.table = FALSE)
    meta <- rbind(meta, meta[1, ])
    meta$sample[nrow(meta)] <- "not_in_the_table.mzML"
    data.table::fwrite(meta, "metadata.csv")

    expect_error(
      suppressMessages(GetMZmineFeature(p$features, "metadata.csv",
                                        group_col = "group", sample_col = "sample")),
      "do not have matching 'area' columns"
    )
    expect_warning(
      suppressMessages(GetMZmineFeature(p$features, "metadata.csv",
                                        group_col = "group", sample_col = "sample",
                                        drop_missing_samples = TRUE)),
      "Proceeding because"
    )
  })
})

test_that("GetMZmineFeature rejects a metadata file without the named columns", {
  p <- pipeline_paths()
  expect_error(
    suppressMessages(GetMZmineFeature(p$features, p$metadata,
                                      group_col = "no_such_column", sample_col = "sample")),
    "group_col must be provided"
  )
  expect_error(
    suppressMessages(GetMZmineFeature(p$features, p$metadata,
                                      group_col = "group", sample_col = "no_such_column")),
    "sample_col must be provided"
  )
})
