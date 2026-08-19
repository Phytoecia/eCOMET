# TRUST: regression -- one test per fixed bug recorded in NEWS.md. The correct
# answer is known because the fix is already in place, so these need no golden
# file. A failure here means the fix was undone.

TOL <- 1e-8

# ---------------------------------------------------------------------------
# NEWS 1: Rao's quadratic entropy normalisation
#
# GetFunctionalHillNumber() divided a features x samples matrix by raoQ, which
# holds one value per sample, so R recycled the divisor down the columns and
# only 4.2% of cells were divided by their own sample's Rao Q. The fix is
# sweep(..., 2, raoQ, "/").
#
# Recycling makes a sample's result depend on where its column sits in the
# matrix. Permuting the columns is therefore a direct probe for it.
# ---------------------------------------------------------------------------

test_that("a sample's functional Hill number does not depend on its column position", {
  fd <- toy_feature_data()
  md <- toy_metadata()
  shuffled <- fd[, c("id", "feature", "s3", "s1", "s2")]
  md_shuffled <- md[c(3, 1, 2), ]

  for (q in c(0, 1, 2)) {
    a <- GetFunctionalHillNumber(fd, md, TOY_SIM, q = q)
    b <- GetFunctionalHillNumber(shuffled, md_shuffled, TOY_SIM, q = q)
    expect_equal(a$hill_number[match(b$sample, a$sample)], b$hill_number,
                 tolerance = TOL, info = paste("q =", q))
  }
})

test_that("the same holds on the fixture, where a single Rao Q would be visible", {
  mmo <- mini_mmo()
  sim <- GetSimMat(mmo, "dreams")
  fd <- mmo$feature_data
  md <- mmo$metadata
  ord <- c(9:1)
  shuffled <- fd[, c("id", "feature", md$sample[ord])]

  a <- GetFunctionalHillNumber(fd, md, sim, q = 1)
  b <- GetFunctionalHillNumber(shuffled, md[ord, ], sim, q = 1)
  expect_equal(a$hill_number[match(b$sample, a$sample)], b$hill_number, tolerance = TOL)

  # Rao Q genuinely varies across these samples, so the test has something to
  # detect: a constant divisor would not reproduce a spread this wide.
  expect_gt(diff(range(a$hill_number)), 1e-6)
})

# ---------------------------------------------------------------------------
# NEWS 2: a sample whose abundances sum to zero returns NA
# ---------------------------------------------------------------------------

test_that("an all-zero sample yields NA, not 1 or Inf", {
  mmo <- toy_mmo_zero_sample()
  fd <- mmo$feature_data
  md <- mmo$metadata
  sim <- GetSimMat(mmo, "toy")

  for (q in c(0, 1, 2)) {
    hill <- GetHillNumbers(fd, md, q = q)
    fh <- GetFunctionalHillNumber(fd, md, sim, q = q)
    expect_true(is.na(hill$hill_number[hill$sample == "s4"]), info = paste("Hill q =", q))
    expect_true(is.na(fh$hill_number[fh$sample == "s4"]), info = paste("functional q =", q))
    # the other samples are unaffected
    expect_false(anyNA(fh$hill_number[fh$sample != "s4"]), info = paste("q =", q))
  }
})

# ---------------------------------------------------------------------------
# NEWS 3: one missing phenotype value no longer zeroes every feature
#
# The fixture metadata carries `mass_missing`, which is `mass` with the value of
# one sample blanked.
# ---------------------------------------------------------------------------

screen_models <- c("pearson", "lm", "spearman")

test_that("a missing phenotype value does not flatten every coefficient", {
  mmo <- mini_mmo()
  expect_true(anyNA(mmo$metadata$mass_missing))   # the fixture must exercise this

  for (model in screen_models) {
    res <- suppressMessages(
      ScreenFeaturePhenotypeCorrelation(mmo, phenotype = "mass_missing",
                                        groups = levels(mmo$metadata$group),
                                        model = model)
    )
    expect_false(all(res$coefficient == 0, na.rm = TRUE), info = model)
    expect_false(all(res$p_value == 1, na.rm = TRUE), info = model)
    expect_true(any(res$p_value < 0.5, na.rm = TRUE), info = model)
  }
})

test_that("dropping the missing sample up front equals removing it from the input", {
  mmo <- mini_mmo()
  drop <- mmo$metadata$sample[is.na(mmo$metadata$mass_missing)]

  complete <- mmo
  complete$metadata <- mmo$metadata[!mmo$metadata$sample %in% drop, ]
  complete$metadata$mass <- complete$metadata$mass_missing
  complete$feature_data <- mmo$feature_data[, setdiff(names(mmo$feature_data), drop)]

  for (model in screen_models) {
    with_na <- suppressMessages(
      ScreenFeaturePhenotypeCorrelation(mmo, phenotype = "mass_missing",
                                        groups = levels(mmo$metadata$group), model = model))
    without <- suppressMessages(
      ScreenFeaturePhenotypeCorrelation(complete, phenotype = "mass",
                                        groups = levels(mmo$metadata$group), model = model))
    expect_equal(with_na$coefficient, without$coefficient, tolerance = TOL, info = model)
    expect_equal(with_na$p_value, without$p_value, tolerance = TOL, info = model)
  }
})

test_that("GetPerformanceFeatureCorrelation and ...Regression survive a missing value", {
  mmo <- mini_mmo()
  ids <- mmo$feature_data$id
  dams <- list(ctrl_vs_sl1.up = ids[1:5])

  corr <- suppressMessages(GetPerformanceFeatureCorrelation(
    mmo, phenotype = "mass_missing", groups = levels(mmo$metadata$group),
    DAM.list = dams, comparisons = character(0)))
  expect_false(all(corr$p_value == 1, na.rm = TRUE))

  regr <- suppressMessages(GetPerformanceFeatureRegression(
    mmo, phenotype = "mass_missing", groups = levels(mmo$metadata$group),
    DAM.list = dams, comparisons = character(0)))
  expect_false(all(regr$p_value == 1, na.rm = TRUE))
})

# ---------------------------------------------------------------------------
# NEWS 4: Spearman uses the exact cor.test p-value, not a normal approximation
# ---------------------------------------------------------------------------

test_that("Spearman p-values equal stats::cor.test(method = 'spearman')", {
  mmo <- mini_mmo()
  res <- suppressMessages(
    ScreenFeaturePhenotypeCorrelation(mmo, phenotype = "mass",
                                      groups = levels(mmo$metadata$group),
                                      model = "spearman"))
  feat <- mmo$feature_data
  y <- mmo$metadata$mass[match(names(feat)[-(1:2)], mmo$metadata$sample)]

  ref <- vapply(seq_len(nrow(feat)), function(i) {
    x <- as.numeric(feat[i, -(1:2)])
    if (stats::sd(x) == 0) return(NA_real_)
    suppressWarnings(stats::cor.test(y, x, method = "spearman")$p.value)
  }, numeric(1))

  keep <- !is.na(ref)
  expect_gt(sum(keep), 20)
  expect_equal(res$p_value[keep], ref[keep], tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# NEWS 5: GetGroupMeans row order and unused factor levels
# ---------------------------------------------------------------------------

test_that("GetGroupMeans returns rows sorted by feature id", {
  gm <- suppressMessages(GetGroupMeans(mini_mmo(), normalization = "None"))
  expect_identical(gm$id, sort(gm$id))
  expect_identical(rownames(gm), as.character(seq_len(nrow(gm))))
})

test_that("GetGroupMeans drops an unused factor level instead of emitting an NA column", {
  mmo <- mini_mmo()
  mmo$metadata$group <- factor(as.character(mmo$metadata$group),
                               levels = c(levels(mmo$metadata$group), "never_measured"))
  gm <- suppressMessages(GetGroupMeans(mmo, normalization = "None"))
  expect_false("never_measured" %in% names(gm))
  expect_false(any(vapply(gm[-1], function(col) all(is.na(col)), logical(1))))
})

# ---------------------------------------------------------------------------
# NEWS 6: use_fastcluster defaults to FALSE and must not change the answer
#
# fastcluster::hclust breaks ties between equal distances differently from
# stats::hclust. On tie-free distances the two must agree exactly.
# ---------------------------------------------------------------------------

test_that("use_fastcluster is FALSE by default", {
  for (fn in c("GetAlphaDiversity", "GetBetaDiversity", "GetFaithPD")) {
    expect_identical(formals(get(fn))$use_fastcluster, FALSE, info = fn)
  }
})

test_that("use_fastcluster = TRUE matches FALSE on tie-free distances", {
  skip_if_not_installed("fastcluster")
  mmo <- suppressMessages(AddCustomSim(mini_mmo(), tie_free_sim(), name = "tiefree"))

  d <- as.dist(1 - tie_free_sim())
  expect_equal(anyDuplicated(as.numeric(d)), 0L)   # the fixture must be tie-free

  a <- suppressMessages(GetBetaDiversity(mmo, method = "Gen.Uni", distance = "tiefree",
                                         normalization = "None", use_fastcluster = FALSE))
  b <- suppressMessages(GetBetaDiversity(mmo, method = "Gen.Uni", distance = "tiefree",
                                         normalization = "None", use_fastcluster = TRUE))
  expect_equal(a, b, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# NEWS 7: GenerateHeatmapInputs stops above 10,000 features
# ---------------------------------------------------------------------------

test_that("GenerateHeatmapInputs refuses to densify more than 10,000 features", {
  n <- 10001L
  ids <- as.character(seq_len(n))
  mmo <- list(
    feature_data = data.frame(id = ids, feature = ids,
                              a = rep(1, n), b = rep(2, n),
                              stringsAsFactors = FALSE),
    metadata = data.frame(sample = c("a", "b"), group = factor(c("g1", "g2")),
                          stringsAsFactors = FALSE)
  )
  class(mmo) <- "mmo"
  # A sparse identity keeps the fixture small; the guard reads nrow() only.
  sim <- Matrix::Diagonal(n)
  dimnames(sim) <- list(ids, ids)
  mmo <- suppressMessages(AddCustomSim(mmo, sim, name = "big"))

  expect_error(
    suppressMessages(GenerateHeatmapInputs(mmo, distance = "big")),
    "requires a dense"
  )
})

# ---------------------------------------------------------------------------
# NEWS 8: similarity slots (.sim) versus the superseded dissimilarity slots
# ---------------------------------------------------------------------------

test_that("GetSimMat names the function to call when only a .dissim slot exists", {
  mmo <- mini_mmo()
  dis <- 1 - GetSimMat(mmo, "dreams")
  mmo$dreams.sim <- NULL
  # "dreams" is a reserved name; the warning is expected here because the test
  # deliberately replaces the .sim slot with a .dissim one of the same name.
  mmo <- suppressWarnings(suppressMessages(AddCustomDist(mmo, dis, name = "dreams")))

  expect_true(!is.null(mmo$dreams.dissim))
  expect_error(GetSimMat(mmo, "dreams"), "AddChemSim")
  expect_error(GetSimMat(mmo, "dreams"), "dreams\\.sim")
})

test_that("AddChemSim and AddCustomSim write the .sim slot the new API reads", {
  mmo <- mini_mmo()
  expect_true(is.null(mmo$dreams.dissim))
  expect_false(is.null(mmo$dreams.sim))
  expect_identical(GetSimMat(mmo, "dreams"), mmo$dreams.sim)

  custom <- suppressMessages(AddCustomSim(mmo, tie_free_sim(), name = "tiefree"))
  expect_false(is.null(custom$tiefree.sim))
  expect_identical(GetSimMat(custom, "tiefree"), custom$tiefree.sim)
})
