# TRUST: external -- each expectation compares eCOMET against an independent
# implementation that eCOMET does not call. GetHillNumbers() implements the
# Hill number formulas directly (R/250813_eCOMET_build_V3.R:4616-4655) and never
# touches vegan, and PairwiseComp() computes a vectorised Welch t-test by hand
# (:2603-2612) rather than calling stats::t.test, so both comparisons are real
# checks rather than restatements.
#
# A failure here means eCOMET and the reference disagree; investigate eCOMET
# first.

TOL <- 1e-8

# ---------------------------------------------------------------------------
# GetHillNumbers vs vegan
# ---------------------------------------------------------------------------

hill_input <- function() {
  mmo <- mini_mmo()
  x <- as.matrix(mmo$feature_data[, -(1:2), drop = FALSE])
  rownames(x) <- mmo$feature_data$id
  t(x)   # vegan wants sites in rows
}

test_that("GetHillNumbers q = 0 matches vegan::specnumber", {
  skip_if_not_installed("vegan")
  mmo <- mini_mmo()
  h <- GetHillNumbers(mmo$feature_data, mmo$metadata, q = 0)
  expect_equal(h$hill_number, unname(vegan::specnumber(hill_input())),
               tolerance = TOL)
})

test_that("GetHillNumbers q = 1 matches exp(vegan Shannon)", {
  skip_if_not_installed("vegan")
  mmo <- mini_mmo()
  h <- GetHillNumbers(mmo$feature_data, mmo$metadata, q = 1)
  expect_equal(h$hill_number,
               unname(exp(vegan::diversity(hill_input(), index = "shannon"))),
               tolerance = TOL)
})

test_that("GetHillNumbers q = 2 matches vegan inverse Simpson", {
  skip_if_not_installed("vegan")
  mmo <- mini_mmo()
  h <- GetHillNumbers(mmo$feature_data, mmo$metadata, q = 2)
  expect_equal(h$hill_number,
               unname(vegan::diversity(hill_input(), index = "invsimpson")),
               tolerance = TOL)
})

# ---------------------------------------------------------------------------
# PairwiseComp vs stats::t.test(var.equal = FALSE) + p.adjust("BH")
#
# The manuscript describes a Welch t-test. PairwiseComp() computes the
# Welch-Satterthwaite degrees of freedom itself; stats::t.test is the reference.
# ---------------------------------------------------------------------------

test_that("PairwiseComp reproduces stats::t.test(var.equal = FALSE)", {
  mmo <- suppressMessages(ReplaceZero(mini_mmo(), method = "one"))
  res <- suppressMessages(
    utils::capture.output(mmo <- PairwiseComp(mmo, "ctrl", "sl1", correction = "BH"))
  )

  feat <- mmo$imputed_feature_data
  meta <- mmo$metadata
  s1 <- meta$sample[meta$group == "ctrl"]
  s2 <- meta$sample[meta$group == "sl1"]

  ref_p <- vapply(seq_len(nrow(feat)), function(i) {
    x <- as.numeric(feat[i, s1]); y <- as.numeric(feat[i, s2])
    if (stats::var(x) == 0 && stats::var(y) == 0) return(NA_real_)
    stats::t.test(x, y, var.equal = FALSE)$p.value
  }, numeric(1))
  ref_fc <- vapply(seq_len(nrow(feat)), function(i) {
    log2(mean(as.numeric(feat[i, s2])) / mean(as.numeric(feat[i, s1])))
  }, numeric(1))

  keep <- !is.na(ref_p)
  expect_gt(sum(keep), 20)   # the comparison must actually cover the fixture

  expect_equal(mmo$pairwise[["ctrl_vs_sl1_pval"]][keep], ref_p[keep],
               tolerance = 1e-10)
  expect_equal(mmo$pairwise[["ctrl_vs_sl1_log2FC"]][keep], ref_fc[keep],
               tolerance = 1e-10)
})

test_that("PairwiseComp adjusts p-values with the requested method", {
  mmo <- suppressMessages(ReplaceZero(mini_mmo(), method = "one"))
  suppressMessages(
    utils::capture.output(mmo <- PairwiseComp(mmo, "ctrl", "sl1", correction = "BH"))
  )
  p <- mmo$pairwise[["ctrl_vs_sl1_pval"]]
  expect_equal(mmo$pairwise[["ctrl_vs_sl1_padj"]],
               stats::p.adjust(p, method = "BH"), tolerance = TOL)
})

test_that("PairwiseComp uses Welch, not Student, degrees of freedom", {
  # Groups with markedly unequal variances separate the two tests. Built by
  # hand so the expected value is unambiguous.
  fd <- data.frame(
    id = "f1", feature = "100.0_1.0",
    a1 = 1, a2 = 2, a3 = 3,
    b1 = 10, b2 = 40, b3 = 100,
    stringsAsFactors = FALSE
  )
  mmo <- list(
    feature_data = fd,
    imputed_feature_data = fd,
    pairwise = data.frame(feature = "100.0_1.0", id = "f1", stringsAsFactors = FALSE),
    metadata = data.frame(
      sample = c("a1", "a2", "a3", "b1", "b2", "b3"),
      group  = factor(rep(c("A", "B"), each = 3)),
      stringsAsFactors = FALSE
    )
  )
  class(mmo) <- "mmo"
  suppressMessages(utils::capture.output(mmo <- PairwiseComp(mmo, "A", "B")))

  x <- c(1, 2, 3); y <- c(10, 40, 100)
  expect_equal(mmo$pairwise[["A_vs_B_pval"]],
               stats::t.test(x, y, var.equal = FALSE)$p.value, tolerance = 1e-10)
  expect_false(isTRUE(all.equal(mmo$pairwise[["A_vs_B_pval"]],
                                stats::t.test(x, y, var.equal = TRUE)$p.value)))
})
