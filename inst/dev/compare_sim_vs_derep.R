# compare_sim_vs_derep.R
#
# Dev-only script — NOT part of the installed package.
#
# PURPOSE:
#   Numerically verify that the new similarity-based functions
#   (using .sim slots / GetSimMat) produce results equivalent to the
#   deprecated _derep functions (using .dissim slots / GetDistanceMat).
#   Uses Tutorial 1 demo data (Arabidopsis treatment study).
#
# FUNCTIONS COMPARED:
#   GetFunctionalHillNumber_derep  vs  GetFunctionalHillNumber
#   GetFaithPD_derep               vs  GetFaithPD
#   GetAlphaDiversity_derep        vs  GetAlphaDiversity   (modes: weighted, faith)
#   GetBetaDiversity_derep         vs  GetBetaDiversity    (methods: CSCS, Gen.Uni, bray, jaccard)
#   GenerateHeatmapInputs_derep    vs  GenerateHeatmapInputs
#
# NOTE ON TREE-BASED FUNCTIONS (GetFaithPD, GetBetaDiversity Gen.Uni):
#   Both _derep and new versions use average-linkage hclust when n <= 10,000,
#   so results SHOULD match on this small demo dataset. If use_mst = TRUE,
#   single-linkage is used and results will intentionally differ.
#
# HOW TO RUN:
#   1. devtools::load_all() from the package root (or install the package)
#   2. source("inst/dev/compare_sim_vs_derep.R")
#   Results print to console. A final summary lists PASS / FAIL for each test.

library(ecomet)
library(dplyr)

# -----------------------------------------------------------------------
# 0. Helper: compare two numeric outputs and report PASS / FAIL
# -----------------------------------------------------------------------
results <- list()

check <- function(label, x, y, tol = 1e-6) {
  # x, y can be numeric vectors, matrices, or data frames (numeric cols only)
  extract_vals <- function(v) {
    if (is.data.frame(v)) {
      num_cols <- sapply(v, is.numeric)
      unlist(v[, num_cols, drop = FALSE])
    } else if (is.list(v)) {
      unlist(lapply(v, function(el) if (is.numeric(el)) el else NULL))
    } else {
      as.numeric(v)
    }
  }
  xv <- extract_vals(x)
  yv <- extract_vals(y)

  eq <- isTRUE(all.equal(xv, yv, tolerance = tol, check.names = FALSE))
  status <- if (eq) "PASS" else "FAIL"
  msg <- if (!eq) paste0("  -> ", all.equal(xv, yv, tolerance = tol, check.names = FALSE)[1]) else ""
  cat(sprintf("[%s]  %s%s\n", status, label, msg))
  results[[label]] <<- eq
  invisible(eq)
}

# -----------------------------------------------------------------------
# 1. Build mmo — Tutorial 1 demo data
# -----------------------------------------------------------------------
cat("\n=== Setup: loading Tutorial 1 demo data ===\n")

data_dir <- system.file("extdata/tutorials/treatment_based", package = "ecomet")
stopifnot(nzchar(data_dir))

mmo <- GetMZmineFeature(
  mzmine_dir   = file.path(data_dir, "feature_table_demo.csv"),
  metadata_dir = file.path(data_dir, "metadata_demo.csv"),
  group_col    = "group",
  sample_col   = "sample"
)
mmo <- LogNormalization(mmo)

# Add BOTH .dissim (for _derep) and .sim (for new functions) to the same mmo
demo_dreams <- file.path(data_dir, "dreams_sim_demo.csv")
mmo <- AddChemDist(mmo, dreams_dir = demo_dreams)  # -> mmo$dreams.dissim
mmo <- AddChemSim(mmo,  dreams_dir = demo_dreams)  # -> mmo$dreams.sim

cat("Features: ", nrow(mmo$feature_data), "\n")
cat("Samples:  ", nrow(mmo$metadata), "\n")
cat("dreams.dissim range: [", min(mmo$dreams.dissim), ",", max(mmo$dreams.dissim), "]\n")
cat("dreams.sim range:    [", min(mmo$dreams.sim),    ",", max(mmo$dreams.sim),    "]\n")
cat("Sum check (dissim + sim diagonal):  ",
    all(abs(diag(mmo$dreams.dissim) + diag(mmo$dreams.sim) - 1) < 1e-10), "\n")

# -----------------------------------------------------------------------
# 2. Extract aligned feature / matrix objects for low-level function tests
# -----------------------------------------------------------------------
feature   <- GetNormFeature(mmo, "Log")
metadata  <- mmo$metadata
dist_mat  <- GetDistanceMat(mmo, "dreams")   # dissimilarity .dissim
sim_mat   <- GetSimMat(mmo,      "dreams")   # similarity    .sim

# Align to features present in both the feature table and the matrix
keep_ids  <- intersect(feature$id, rownames(dist_mat))
feat_sub  <- feature[feature$id %in% keep_ids, , drop = FALSE]
dist_sub  <- dist_mat[keep_ids, keep_ids]
sim_sub   <- sim_mat[keep_ids, keep_ids]

cat("\nAligned feature count: ", length(keep_ids), "\n")
cat("D + S == 1 check (spot):  ",
    isTRUE(all.equal(as.numeric(dist_sub[1:5, 1:5] + sim_sub[1:5, 1:5]),
                     rep(1, 25), tolerance = 1e-6)), "\n")

# -----------------------------------------------------------------------
# 3. GetFunctionalHillNumber_derep  vs  GetFunctionalHillNumber
# -----------------------------------------------------------------------
cat("\n=== GetFunctionalHillNumber ===\n")

hill_derep_q1 <- GetFunctionalHillNumber_derep(
  feature = feat_sub, metadata = metadata, distance_matrix = dist_sub, q = 1
)
hill_new_q1 <- GetFunctionalHillNumber(
  feature = feat_sub, metadata = metadata, sim_matrix = sim_sub, q = 1
)
check("FunctionalHillNumber q=1  (value column)", hill_derep_q1$value, hill_new_q1$value)

hill_derep_q2 <- GetFunctionalHillNumber_derep(
  feature = feat_sub, metadata = metadata, distance_matrix = dist_sub, q = 2
)
hill_new_q2 <- GetFunctionalHillNumber(
  feature = feat_sub, metadata = metadata, sim_matrix = sim_sub, q = 2
)
check("FunctionalHillNumber q=2  (value column)", hill_derep_q2$value, hill_new_q2$value)

# -----------------------------------------------------------------------
# 4. GetFaithPD_derep  vs  GetFaithPD  (average linkage, n <= 10k)
# -----------------------------------------------------------------------
cat("\n=== GetFaithPD (average linkage — n <= 10,000) ===\n")

fpd_derep <- GetFaithPD_derep(
  feature = feat_sub, metadata = metadata, distance_matrix = dist_sub
)
fpd_new <- GetFaithPD(
  feature = feat_sub, metadata = metadata, sim_matrix = sim_sub,
  use_mst = FALSE   # use average linkage to match _derep
)
check("FaithPD               (PD column)",  fpd_derep$PD, fpd_new$PD)
check("FaithPD               (SR column)",  fpd_derep$SR, fpd_new$SR)

cat("\n--- FaithPD with use_mst = TRUE (single linkage — different tree, values WILL differ) ---\n")
fpd_mst <- GetFaithPD(
  feature = feat_sub, metadata = metadata, sim_matrix = sim_sub,
  use_mst = TRUE
)
cat("  PD range _derep: [", round(min(fpd_derep$PD), 4), ",", round(max(fpd_derep$PD), 4), "]\n")
cat("  PD range MST:    [", round(min(fpd_mst$PD),   4), ",", round(max(fpd_mst$PD),   4), "]\n")
cat("  (Difference is expected — single linkage builds a different tree)\n")

# -----------------------------------------------------------------------
# 5. GetAlphaDiversity_derep  vs  GetAlphaDiversity
# -----------------------------------------------------------------------
cat("\n=== GetAlphaDiversity ===\n")

for (md in c("weighted", "faith")) {
  for (q_val in c(1, 2)) {
    if (md == "faith" && q_val == 2) next  # faith has no q argument

    alpha_derep <- GetAlphaDiversity_derep(
      mmo, mode = md, q = q_val, normalization = "Log",
      distance = "dreams", output = "sample_level"
    )
    alpha_new <- GetAlphaDiversity(
      mmo, mode = md, q = q_val, normalization = "Log",
      distance = "dreams", output = "sample_level", use_mst = FALSE
    )
    lbl <- sprintf("AlphaDiversity mode=%-8s q=%d  (value column)", md, q_val)
    check(lbl, alpha_derep$value, alpha_new$value)
  }
}

# -----------------------------------------------------------------------
# 6. GetBetaDiversity_derep  vs  GetBetaDiversity
# -----------------------------------------------------------------------
cat("\n=== GetBetaDiversity ===\n")

# CSCS: algebraically equivalent (CSS = S directly vs CSS = 1 - D)
beta_cscs_derep <- GetBetaDiversity_derep(
  mmo, method = "CSCS", normalization = "Log", distance = "dreams"
)
beta_cscs_new <- GetBetaDiversity(
  mmo, method = "CSCS", normalization = "Log", distance = "dreams"
)
check("BetaDiversity CSCS    (full dissimilarity matrix)", beta_cscs_derep, beta_cscs_new)

# Gen.Uni: both average-linkage at small n — should match
beta_genuni_derep <- GetBetaDiversity_derep(
  mmo, method = "Gen.Uni", normalization = "Log", distance = "dreams"
)
beta_genuni_new <- GetBetaDiversity(
  mmo, method = "Gen.Uni", normalization = "Log", distance = "dreams",
  use_mst = FALSE
)
# Gen.Uni returns a list (d_0, d_0.5, d_1); compare each
for (nm in names(beta_genuni_derep)) {
  if (is.matrix(beta_genuni_derep[[nm]])) {
    check(
      sprintf("BetaDiversity Gen.Uni (%s)", nm),
      beta_genuni_derep[[nm]], beta_genuni_new[[nm]]
    )
  }
}

# Bray-Curtis and Jaccard: no distance matrix at all — must be identical
beta_bray_derep <- GetBetaDiversity_derep(mmo, method = "bray",    normalization = "Log")
beta_bray_new   <- GetBetaDiversity(mmo,       method = "bray",    normalization = "Log")
check("BetaDiversity bray     (full matrix)", beta_bray_derep, beta_bray_new)

beta_jacc_derep <- GetBetaDiversity_derep(mmo, method = "jaccard", normalization = "PA")
beta_jacc_new   <- GetBetaDiversity(mmo,       method = "jaccard", normalization = "PA")
check("BetaDiversity jaccard  (full matrix)", beta_jacc_derep, beta_jacc_new)

# -----------------------------------------------------------------------
# 7. GenerateHeatmapInputs_derep  vs  GenerateHeatmapInputs
# -----------------------------------------------------------------------
cat("\n=== GenerateHeatmapInputs ===\n")

hm_derep <- GenerateHeatmapInputs_derep(
  mmo, normalization = "Log", distance = "dreams"
)
hm_new <- GenerateHeatmapInputs(
  mmo, normalization = "Log", distance = "dreams"
)
# Both return a list; compare the main data matrix (heatmap_data)
check("HeatmapInputs         (heatmap_data matrix)",
      hm_derep$heatmap_data, hm_new$heatmap_data)
# Clustering row/col order may differ by floating-point ties; compare the
# distance matrix used for clustering instead
if (!is.null(hm_derep$dist_matrix) && !is.null(hm_new$dist_matrix)) {
  check("HeatmapInputs         (dist_matrix used for clustering)",
        hm_derep$dist_matrix, hm_new$dist_matrix)
}

# -----------------------------------------------------------------------
# 8. Summary
# -----------------------------------------------------------------------
cat("\n", strrep("=", 60), "\n", sep = "")
cat("SUMMARY\n")
cat(strrep("=", 60), "\n", sep = "")
pass <- sum(unlist(results))
fail <- sum(!unlist(results))
cat(sprintf("PASS: %d / %d\n", pass, pass + fail))
if (fail > 0) {
  cat("FAIL:\n")
  for (nm in names(results)) if (!results[[nm]]) cat("  -", nm, "\n")
}
cat(strrep("=", 60), "\n", sep = "")
