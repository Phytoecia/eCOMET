# helper-mmo.R --------------------------------------------------------------
#
# Fixture loaders and the hand-checkable toy object used across the suite.
#
# The fixtures under fixtures/ are produced by data-raw/make_test_fixture.R
# from the bundled tutorial data. Regenerate them with
#   Rscript data-raw/make_test_fixture.R
# after any change to the importers; never edit them by hand.
# ---------------------------------------------------------------------------

# Absolute, so a path captured before withr::with_tempdir() still resolves
# after the working directory changes.
fixture_path <- function(...) {
  normalizePath(testthat::test_path("fixtures", ...), mustWork = TRUE)
}

# Loaders are memoised so a full test run deserialises each fixture once.
.fixture_cache <- new.env(parent = emptyenv())

.load_fixture <- function(file) {
  if (is.null(.fixture_cache[[file]])) {
    .fixture_cache[[file]] <- readRDS(fixture_path(file))
  }
  .fixture_cache[[file]]
}

#' 30 annotated features x 9 samples (3 groups of 3), with a DreaMS similarity
#' matrix in mmo$dreams.sim. Every feature carries a CANOPUS annotation.
mini_mmo <- function() .load_fixture("mini_mmo.rds")

#' The same slice plus three features that SIRIUS did not annotate, for the
#' missing-annotation paths.
mini_mmo_unannotated <- function() .load_fixture("mini_mmo_with_unannotated.rds")


# -- toy object -------------------------------------------------------------
#
# Three features x three samples with abundances chosen so that every quantity
# the `hand` tests assert has a closed form.
#
#   proportions      p1 = (1/3, 1/3, 1/3)
#                    p2 = (1/2, 1/4,  1/4)
#                    p3 = (3/5, 1/5,  1/5)
#
#   similarity S     f1-f2 = 1/2,  f1-f3 = 0,  f2-f3 = 1/4,  diagonal 1
#
#   Rao Q = 1 - p'Sp  ->  Q1 = 1/2,  Q2 = 0.53125,  Q3 = 0.42
#
# See test-diversity-hand.R for the derivations that use these numbers.

TOY_SIM <- local({
  s <- matrix(c(1,   0.5, 0,
                0.5, 1,   0.25,
                0,   0.25, 1), nrow = 3, byrow = TRUE)
  dimnames(s) <- list(c("f1", "f2", "f3"), c("f1", "f2", "f3"))
  s
})

toy_feature_data <- function() {
  data.frame(
    id      = c("f1", "f2", "f3"),
    feature = c("100.0_1.0", "200.0_2.0", "300.0_3.0"),
    s1 = c(10, 10, 10),
    s2 = c(20, 10, 10),
    s3 = c(60, 20, 20),
    stringsAsFactors = FALSE
  )
}

toy_metadata <- function() {
  data.frame(
    sample = c("s1", "s2", "s3"),
    group  = factor(c("A", "A", "B"), levels = c("A", "B")),
    pheno  = c(1.5, 2.5, 4.0),
    stringsAsFactors = FALSE
  )
}

#' A minimal but complete mmo built from the toy numbers above.
#'
#' The similarity matrix is stored as `toy.sim`, so downstream calls use
#' `distance = "toy"`. AddCustomSim() warns for the reserved built-in names.
toy_mmo <- function() {
  fd <- toy_feature_data()
  mmo <- list(
    feature_data = fd,
    feature_info = data.frame(
      id = fd$id, feature = fd$feature,
      mz = c(100, 200, 300), rt = c(1, 2, 3),
      stringsAsFactors = FALSE
    ),
    pairwise = data.frame(feature = fd$feature, id = fd$id, stringsAsFactors = FALSE),
    metadata = toy_metadata()
  )
  class(mmo) <- "mmo"
  mmo <- FeaturePresence(mmo, threshold = 1)
  suppressMessages(AddCustomSim(mmo, sim_matrix = TOY_SIM, name = "toy"))
}

#' The toy object with a fourth sample whose abundances are all zero.
#'
#' Used for the NEWS regression that an empty sample yields NA rather than 1
#' (q = 1) or Inf (q != 1).
toy_mmo_zero_sample <- function() {
  mmo <- toy_mmo()
  mmo$feature_data$s4 <- c(0, 0, 0)
  mmo$metadata <- rbind(
    mmo$metadata,
    data.frame(sample = "s4", group = factor("B", levels = c("A", "B")),
               pheno = 3.0, stringsAsFactors = FALSE)
  )
  mmo <- FeaturePresence(mmo, threshold = 1)
  mmo
}

#' A tie-free similarity matrix over the fixture's features.
#'
#' The DreaMS scores in the fixture are rounded to three decimals, so the
#' derived distances contain ties and hierarchical clustering is only
#' deterministic within one implementation. Tests that compare `stats::hclust`
#' against `fastcluster::hclust` need distances that are all distinct.
tie_free_sim <- function(mmo = mini_mmo()) {
  ids <- mmo$feature_data$id
  n <- length(ids)
  # Distinct, reproducible off-diagonal values in (0, 1): no RNG involved.
  vals <- seq_len(n * (n - 1) / 2) / (n * (n - 1) / 2 + 1)
  s <- matrix(0, n, n, dimnames = list(ids, ids))
  s[lower.tri(s)] <- vals
  s <- s + t(s)
  diag(s) <- 1
  s
}
