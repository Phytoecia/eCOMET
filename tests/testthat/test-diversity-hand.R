# TRUST: hand -- every expectation in this file is a value derived on paper from
# the toy numbers in helper-mmo.R, not a recording of what the code currently
# returns. A failure here means the arithmetic no longer matches the published
# definition, i.e. a bug in eCOMET, and must not be "fixed" by editing the
# expected value.

TOL <- 1e-8

# ---------------------------------------------------------------------------
# GetHillNumbers
#
#   p1 = (1/3, 1/3, 1/3)   p2 = (1/2, 1/4, 1/4)   p3 = (3/5, 1/5, 1/5)
#
#   q = 0  richness                 -> 3 for every sample
#   q = 1  exp(Shannon)             -> 3, 2^(3/2), 1 / (0.6^0.6 * 0.2^0.4)
#   q = 2  inverse Simpson 1/sum p^2-> 3, 8/3, 25/11
# ---------------------------------------------------------------------------

test_that("GetHillNumbers q = 0 is observed richness", {
  h <- GetHillNumbers(toy_feature_data(), toy_metadata(), q = 0)
  expect_equal(h$hill_number, c(3, 3, 3), tolerance = TOL)
})

test_that("GetHillNumbers q = 1 is exp(Shannon)", {
  h <- GetHillNumbers(toy_feature_data(), toy_metadata(), q = 1)
  expected <- c(
    3,                                   # uniform over three features
    2^(3/2),                             # -(1/2 log 1/2 + 2 * 1/4 log 1/4) = (3/2) log 2
    1 / (0.6^0.6 * 0.2^0.2 * 0.2^0.2)    # exp(-sum p log p) written as a product
  )
  expect_equal(h$hill_number, expected, tolerance = TOL)
})

test_that("GetHillNumbers q = 2 is inverse Simpson", {
  h <- GetHillNumbers(toy_feature_data(), toy_metadata(), q = 2)
  expected <- c(
    3,          # 1 / (3 * (1/3)^2)
    8 / 3,      # 1 / (1/4 + 1/16 + 1/16)
    25 / 11     # 1 / (0.36 + 0.04 + 0.04)
  )
  expect_equal(h$hill_number, expected, tolerance = TOL)
})

test_that("GetHillNumbers reports the sample and group of each column", {
  h <- GetHillNumbers(toy_feature_data(), toy_metadata(), q = 1)
  expect_equal(h$sample, c("s1", "s2", "s3"))
  expect_equal(as.character(h$group), c("A", "A", "B"))
})

# ---------------------------------------------------------------------------
# GetFunctionalHillNumber
#
# Rao Q = 1 - p'Sp:
#   Q1 = 1 - 4.5/9      = 1/2
#   Q2 = 1 - 0.53125    = 0.46875
#   Q3 = 1 - 0.58       = 0.42
#
# For an assemblage with uniform proportions every order of the functional Hill
# number collapses to the same value, Q / p^2 (all Hill numbers coincide when
# the distribution is uniform). With p = 1/3 and Q1 = 1/2 that is
#
#   0.5 * 9 = 4.5    at q = 0, 1 and 2.
# ---------------------------------------------------------------------------

test_that("GetFunctionalHillNumber of a uniform sample equals Q / p^2 at every q", {
  fd <- toy_feature_data()[, c("id", "feature", "s1")]
  md <- toy_metadata()[1, ]
  for (q in c(0, 1, 2)) {
    fh <- GetFunctionalHillNumber(fd, md, TOY_SIM, q = q)
    expect_equal(fh$hill_number, 4.5, tolerance = TOL,
                 info = paste("q =", q))
  }
})

test_that("each sample is divided by its own Rao Q (NEWS: sweep, not recycling)", {
  # The recycling bug divided the whole feature x sample proportion matrix by a
  # per-sample vector, so a feature's proportion could be divided by another
  # sample's Rao Q. Computing one sample at a time removes any opportunity to
  # recycle: with the sweep() fix the two routes must agree exactly.
  fd <- toy_feature_data()
  md <- toy_metadata()
  for (q in c(0, 1, 2)) {
    together <- GetFunctionalHillNumber(fd, md, TOY_SIM, q = q)
    alone <- vapply(md$sample, function(s) {
      GetFunctionalHillNumber(fd[, c("id", "feature", s)],
                              md[md$sample == s, ], TOY_SIM, q = q)$hill_number
    }, numeric(1))
    expect_equal(together$hill_number, unname(alone), tolerance = TOL,
                 info = paste("q =", q))
  }
})

test_that("the Rao Q implied by GetFunctionalHillNumber matches the hand value", {
  # At q = 2 the implemented formula reduces to
  #   FHN = ( sum_i p_i^2 (sum_j p_j^2 - (S p^2)_i) / Q^2 )^-1
  # so a single feature-uniform sample lets Q be read back off the result:
  # FHN = Q / p^2  =>  Q = FHN * p^2.
  fh <- GetFunctionalHillNumber(toy_feature_data()[, c("id", "feature", "s1")],
                                toy_metadata()[1, ], TOY_SIM, q = 2)
  expect_equal(fh$hill_number * (1/3)^2, 0.5, tolerance = TOL)
})

# ---------------------------------------------------------------------------
# GetBetaDiversity(method = "CSCS")
#
#   C = P' S P with columns of P the sample proportions:
#     C11 = 1/2        C12 = 1/2       C13 = 1/2
#                      C22 = 0.53125   C23 = 0.55
#                                      C33 = 0.58
#   CSCS[i,j] = C[i,j] / max(C[i,i], C[j,j]);  beta = 1 - CSCS.
# ---------------------------------------------------------------------------

test_that("CSCS beta diversity matches the hand-computed matrix", {
  beta <- suppressMessages(
    GetBetaDiversity(toy_mmo(), method = "CSCS", distance = "toy",
                     normalization = "None")
  )
  expected <- matrix(0, 3, 3, dimnames = list(c("s1", "s2", "s3"),
                                              c("s1", "s2", "s3")))
  expected["s1", "s2"] <- expected["s2", "s1"] <- 1 - 0.5  / 0.53125
  expected["s1", "s3"] <- expected["s3", "s1"] <- 1 - 0.5  / 0.58
  expected["s2", "s3"] <- expected["s3", "s2"] <- 1 - 0.55 / 0.58
  expect_equal(beta, expected, tolerance = TOL)
})

test_that("CSCS is zero on the diagonal and symmetric", {
  beta <- suppressMessages(
    GetBetaDiversity(toy_mmo(), method = "CSCS", distance = "toy",
                     normalization = "None")
  )
  expect_equal(diag(beta), c(s1 = 0, s2 = 0, s3 = 0), tolerance = TOL)
  expect_true(isSymmetric(beta))
})

# ---------------------------------------------------------------------------
# CanopusLevelEnrichmentAnal
#
# Ten annotated features: pathway A x 4, B x 3, C x 3. The tested list holds
# three A features and one B feature. CanopusLevelEnrichmentAnal builds
#
#         term   not-term
#   list    3        1        (row total 4  = size of the tested list)
#   all     4        6        (row total 10 = all annotated features)
#
# and passes it to fisher.test(alternative = "greater"). Column totals are 7
# and 7, grand total 14, so under the hypergeometric null
#
#   P(X >= 3) = [C(7,3) C(7,1) + C(7,4) C(7,0)] / C(14,4)
#             = (245 + 35) / 1001 = 280 / 1001.
#
# Fold enrichment = (3/4) / (4/10) = 1.875.
#
# Note that the second row counts ALL features, not the features outside the
# tested list, so the two rows overlap. The expectations below record the table
# the code actually builds.
# ---------------------------------------------------------------------------

enrichment_toy_mmo <- function() {
  ids <- paste0("f", 1:10)
  mmo <- list(
    feature_data = data.frame(id = ids, feature = ids, stringsAsFactors = FALSE),
    sirius_annot = data.frame(
      id = ids,
      `NPC#pathway` = c(rep("A", 4), rep("B", 3), rep("C", 3)),
      check.names = FALSE, stringsAsFactors = FALSE
    )
  )
  class(mmo) <- "mmo"
  mmo
}

test_that("CanopusLevelEnrichmentAnal reproduces the hand-built 2x2 table", {
  res <- CanopusLevelEnrichmentAnal(
    enrichment_toy_mmo(), list_test = c("f1", "f2", "f3", "f5"),
    term_level = "NPC_pathway", representation = "greater", sig = FALSE
  )
  res <- res[order(res$term), ]

  expect_equal(res$term, c("A", "B"))
  expect_equal(res$subsetcount, c(3, 1))
  expect_equal(res$totalcount, c(4, 3))
  expect_equal(res$foldenrichment, c((3/4) / (4/10), (1/4) / (3/10)), tolerance = TOL)
  expect_equal(res$pval[res$term == "A"], 280 / 1001, tolerance = TOL)
})

test_that("CanopusLevelEnrichmentAnal rejects an unknown term level", {
  expect_error(
    CanopusLevelEnrichmentAnal(enrichment_toy_mmo(), list_test = "f1",
                               term_level = "not_a_level"),
    "Invalid term level"
  )
})
