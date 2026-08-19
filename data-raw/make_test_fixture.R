# make_test_fixture.R -------------------------------------------------------
#
# Builds the small fixtures used by tests/testthat from the bundled tutorial
# data. Run manually from the package root:
#
#   Rscript data-raw/make_test_fixture.R
#
# The tutorial data under inst/extdata/tutorials/ is 82 MB; tests must never
# read it. This script cuts a ~30-feature x 9-sample slice that preserves the
# properties the tests depend on:
#
#   * every fixture feature carries a SIRIUS/CANOPUS annotation, because the
#     annotation rate of a random slice would leave the annotation-dependent
#     functions untested;
#   * at least three NPC pathways are represented, so the class-enrichment
#     contingency tables are non-degenerate;
#   * some features hold zeros (ReplaceZero) and some samples hold missing
#     phenotype values (the NEWS missing-phenotype regression);
#   * mini_feature_table.csv keeps the MZmine column headers verbatim, so the
#     integration test breaks if an MZmine export format changes.
#
# This script only READS the tutorial data. It never writes into inst/extdata.
#
# Revision history:
#   2026-08-19  initial version (test suite for NPH revision)
# ---------------------------------------------------------------------------

suppressMessages({
  library(data.table)
  devtools::load_all(".", quiet = TRUE)
})

set.seed(20260819)

SRC     <- "inst/extdata/tutorials/treatment_based"
FIXTURE <- "tests/testthat/fixtures"
dir.create(FIXTURE, recursive = TRUE, showWarnings = FALSE)

N_FEATURES     <- 30L   # annotated features in the main fixture
N_UNANNOTATED   <- 3L   # extra unannotated features, kept out of the main fixture
SAMPLES <- c("con_1", "con_2", "con_3",
             "sl1_1", "sl1_2", "sl1_3",
             "le1_1", "le1_2", "le1_3")

msg <- function(...) cat("[fixture]", ..., "\n")

# -- 1. read the tutorial sources ------------------------------------------
ft <- fread(file.path(SRC, "feature_table_demo.csv"), check.names = FALSE,
            na.strings = c("", "NA"), data.table = FALSE)
can <- fread(file.path(SRC, "canopus_formula_summary.tsv"), check.names = FALSE,
             na.strings = c("", "NA"), data.table = FALSE)
str_id <- fread(file.path(SRC, "structure_identifications.tsv"), check.names = FALSE,
                na.strings = c("", "NA"), data.table = FALSE)
meta <- fread(file.path(SRC, "metadata_demo.csv"), check.names = FALSE,
              na.strings = c("", "NA"), data.table = FALSE)

msg("sources:", nrow(ft), "features,", nrow(can), "CANOPUS rows,",
    nrow(str_id), "structure rows,", nrow(meta), "samples")

ft_id  <- gsub(" ", "", as.character(ft$id))
can_id <- gsub(" ", "", as.character(can$mappingFeatureId))
str_ids <- gsub(" ", "", as.character(str_id$mappingFeatureId))

# -- 2. choose features -----------------------------------------------------
# Area columns for the nine retained samples, in SAMPLES order.
area_cols <- paste0("datafile:", SAMPLES, ".mzML:area")
stopifnot(all(area_cols %in% names(ft)))

area <- as.matrix(ft[, area_cols, drop = FALSE])
area[is.na(area)] <- 0
rownames(area) <- ft_id
colnames(area) <- SAMPLES

n_zero  <- rowSums(area == 0)
total   <- rowSums(area)

annotated <- ft_id %in% can_id
# Usable candidates: annotated, detected in at least five of the nine samples,
# so that relative proportions and Rao Q stay well defined.
cand <- annotated & n_zero <= 4L & total > 0

pathway <- can[["NPC#pathway"]][match(ft_id, can_id)]
npc_class <- can[["NPC#class"]][match(ft_id, can_id)]

pool <- data.frame(
  id       = ft_id,
  pathway  = pathway,
  npc_class = npc_class,
  n_zero   = n_zero,
  total    = total,
  in_structure = ft_id %in% str_ids,
  stringsAsFactors = FALSE
)[cand, ]

# Deterministic order: by pathway, structure-annotated first, then numeric id.
# Preferring features that also appear in structure_identifications.tsv keeps
# the COSMIC/structure columns of mmo$sirius_annot populated for most rows.
pool <- pool[order(pool$pathway, !pool$in_structure, as.numeric(pool$id)), ]

# Take features with zeros first (ReplaceZero needs them), then fill the rest
# round-robin across pathways so no single NPC pathway dominates.
pick_round_robin <- function(df, n) {
  if (n <= 0 || nrow(df) == 0) return(character(0))
  by_pw <- split(df$id, df$pathway)
  out <- character(0)
  k <- 1L
  while (length(out) < n && any(lengths(by_pw) > 0)) {
    for (pw in names(by_pw)) {
      if (length(out) >= n) break
      if (length(by_pw[[pw]]) >= k) out <- c(out, by_pw[[pw]][k])
    }
    k <- k + 1L
  }
  head(out, n)
}

with_zero <- pool[pool$n_zero > 0, ]
no_zero   <- pool[pool$n_zero == 0, ]

N_WITH_ZERO <- 6L
sel_zero <- pick_round_robin(with_zero, N_WITH_ZERO)
sel_rest <- pick_round_robin(no_zero[!no_zero$id %in% sel_zero, ],
                             N_FEATURES - length(sel_zero))
sel_annot <- c(sel_zero, sel_rest)
stopifnot(length(sel_annot) == N_FEATURES, !anyDuplicated(sel_annot))

# Unannotated extras for the missing-annotation fixture.
unannot_pool <- ft_id[!annotated & n_zero <= 4L & total > 0]
sel_unannot <- head(unannot_pool[order(as.numeric(unannot_pool))], N_UNANNOTATED)
stopifnot(length(sel_unannot) == N_UNANNOTATED)

sel_all <- c(sel_annot, sel_unannot)

msg("selected", length(sel_annot), "annotated +", length(sel_unannot), "unannotated features")
msg("NPC pathways covered:", length(unique(pool$pathway[pool$id %in% sel_annot])))
print(table(pool$pathway[pool$id %in% sel_annot]))
msg("NPC classes covered:", length(unique(pool$npc_class[pool$id %in% sel_annot])))
msg("features containing a zero:", sum(pool$n_zero[pool$id %in% sel_annot] > 0))
msg("features also in structure_identifications:",
    sum(pool$in_structure[pool$id %in% sel_annot]))

# -- 3. mini_feature_table.csv ---------------------------------------------
# Keep every MZmine feature-level column and the per-sample datafile columns of
# the nine retained samples. Column NAMES are never rewritten: the integration
# test exists to notice when MZmine changes them.
datafile_cols <- grep("^datafile:", names(ft), value = TRUE)
keep_datafile <- grep(paste0("^datafile:(", paste(SAMPLES, collapse = "|"), ")\\.mzML:"),
                      datafile_cols, value = TRUE)
keep_cols <- c(setdiff(names(ft), datafile_cols), keep_datafile)

mini_ft <- ft[match(sel_all, ft_id), keep_cols, drop = FALSE]
fwrite(mini_ft, file.path(FIXTURE, "mini_feature_table.csv"))
msg("mini_feature_table.csv:", nrow(mini_ft), "rows x", ncol(mini_ft), "cols")

# -- 4. mini_metadata.csv ---------------------------------------------------
meta_sample <- names(meta)[1]           # the source file starts with a UTF-8 BOM
mini_meta <- meta[sub("\\.mzML$", "", meta[[meta_sample]]) %in% SAMPLES, , drop = FALSE]
mini_meta <- mini_meta[match(SAMPLES, sub("\\.mzML$", "", mini_meta[[meta_sample]])), ]
# `mass` is complete; `sl` and `le_rgr` are missing outside their own treatment.
# `mass_missing` is `mass` with one value blanked so that a phenotype which is
# mostly present but not complete exists - the shape that triggered the NEWS
# missing-phenotype bug.
mini_meta$mass_missing <- mini_meta$mass
mini_meta$mass_missing[2] <- NA_real_
names(mini_meta)[1] <- "sample"
# The tutorial metadata starts with a UTF-8 BOM. data.table::fread() strips it,
# so GetMZmineFeature() reads the first column as "sample" - keep the BOM here
# so the fixture still exercises that path.
.meta_tmp <- tempfile(fileext = ".csv")
fwrite(mini_meta, .meta_tmp)
.con <- file(file.path(FIXTURE, "mini_metadata.csv"), "wb")
writeBin(as.raw(c(0xEF, 0xBB, 0xBF)), .con)
writeBin(readBin(.meta_tmp, "raw", file.size(.meta_tmp)), .con)
close(.con)
unlink(.meta_tmp)
msg("mini_metadata.csv:", nrow(mini_meta), "samples")
print(mini_meta)

# -- 5. mini_canopus.tsv / mini_structure.tsv ------------------------------
mini_can <- can[can_id %in% sel_all, , drop = FALSE]
fwrite(mini_can, file.path(FIXTURE, "mini_canopus.tsv"), sep = "\t")
mini_str <- str_id[str_ids %in% sel_all, , drop = FALSE]
fwrite(mini_str, file.path(FIXTURE, "mini_structure.tsv"), sep = "\t")
msg("mini_canopus.tsv:", nrow(mini_can), "rows; mini_structure.tsv:", nrow(mini_str), "rows")

# -- 6. mini_dreams_sim.csv -------------------------------------------------
# The DreaMS export is an edge list (ID1, ID2, EdgeType, Score, Annotation).
# Keep the header verbatim and the edges whose both endpoints survive.
sim_edges <- fread(file.path(SRC, "dreams_sim_demo.csv"), check.names = FALSE,
                   data.table = FALSE, colClasses = "character")
hdr <- names(sim_edges)
keep_edge <- trimws(sim_edges[[1]]) %in% sel_all & trimws(sim_edges[[2]]) %in% sel_all
mini_sim <- sim_edges[keep_edge, , drop = FALSE]
fwrite(mini_sim, file.path(FIXTURE, "mini_dreams_sim.csv"))
msg("mini_dreams_sim.csv:", nrow(mini_sim), "edges of",
    length(sel_all) * (length(sel_all) - 1) / 2, "possible")

# -- 6b. mini_custom_db.csv -------------------------------------------------
# A custom compound database in the tutorial's format (compound, mz, rt). The
# entries are the exact mz/rt of five fixture features, so AddCustomAnnot()
# has guaranteed hits inside its default 5 ppm / 0.5 min tolerances.
db_src <- mini_ft[match(head(sel_annot, 5), gsub(" ", "", as.character(mini_ft$id))), ]
mini_db <- data.frame(
  compound = paste0("test_compound_", seq_len(nrow(db_src))),
  mz = db_src$mz,
  rt = db_src$rt,
  stringsAsFactors = FALSE
)
fwrite(mini_db, file.path(FIXTURE, "mini_custom_db.csv"))
msg("mini_custom_db.csv:", nrow(mini_db), "compounds")

# -- 6c. mini_spectra.mgf ---------------------------------------------------
# A SYNTHETIC MGF. The tutorial data ships no MGF, and the two MGF readers only
# parse FEATURE_ID= and MSLEVEL=, so a minimal well-formed file exercises them.
# Only the first 20 features get a spectrum, leaving features without MS2 so
# filter_mgf_to_mmo() and annotate_feature_info_ms2_from_mgf() have something to
# drop and something to report as missing. Peak lists are placeholders.
mgf_ids <- head(sel_all, 20)
mgf_lines <- unlist(lapply(mgf_ids, function(fid) {
  row <- mini_ft[gsub(" ", "", as.character(mini_ft$id)) == fid, ][1, ]
  c("BEGIN IONS",
    paste0("FEATURE_ID=", fid),
    paste0("SCANS=", fid),
    paste0("PEPMASS=", format(row$mz, digits = 10)),
    paste0("RTINSECONDS=", format(round(row$rt * 60, 3), nsmall = 3)),
    "CHARGE=1+",
    "MSLEVEL=2",
    sprintf("%.4f %.1f", row$mz * 0.5, 120.0),
    sprintf("%.4f %.1f", row$mz * 0.75, 540.0),
    sprintf("%.4f %.1f", row$mz, 999.0),
    "END IONS",
    "")
}))
writeLines(mgf_lines, file.path(FIXTURE, "mini_spectra.mgf"))
msg("mini_spectra.mgf:", length(mgf_ids), "spectra")

# -- 7. build the mmo objects ----------------------------------------------
mmo_full <- GetMZmineFeature(
  mzmine_dir   = file.path(FIXTURE, "mini_feature_table.csv"),
  metadata_dir = file.path(FIXTURE, "mini_metadata.csv"),
  group_col    = "group",
  sample_col   = "sample"
)
mmo_full <- AddSiriusAnnot(
  mmo_full,
  canopus_structuredir = file.path(FIXTURE, "mini_structure.tsv"),
  canopus_formuladir   = file.path(FIXTURE, "mini_canopus.tsv")
)
mmo_full <- AddChemSim(mmo_full, dreams_dir = file.path(FIXTURE, "mini_dreams_sim.csv"))

mmo <- filter_mmo(mmo_full, id_list = sel_annot, drop_empty_feat = FALSE)

saveRDS(mmo, file.path(FIXTURE, "mini_mmo.rds"), compress = "xz")
saveRDS(mmo_full, file.path(FIXTURE, "mini_mmo_with_unannotated.rds"), compress = "xz")

# -- 8. QC ------------------------------------------------------------------
msg("--- QC ---")
stopifnot(nrow(mmo$feature_data) == N_FEATURES)
stopifnot(ncol(mmo$feature_data) == length(SAMPLES) + 2L)
stopifnot(nrow(mmo$metadata) == length(SAMPLES))
stopifnot(all(!is.na(mmo$sirius_annot[["NPC#pathway"]])))
stopifnot(length(unique(mmo$sirius_annot[["NPC#pathway"]])) >= 3L)
stopifnot(any(as.matrix(mmo$feature_data[, -(1:2)]) == 0))
stopifnot(anyNA(mmo$metadata$mass_missing))

sim <- GetSimMat(mmo, "dreams")
stopifnot(isSymmetric(unname(as.matrix(sim))))
stopifnot(all(abs(diag(as.matrix(sim)) - 1) < 1e-6))

# DreaMS scores are reported to three decimals, so a 30-feature slice always
# carries tied distances. That is faithful to the real data, and stats::hclust()
# resolves ties deterministically, so golden snapshots stay stable. The
# use_fastcluster equivalence test builds its own tie-free matrix instead of
# relying on this one.
d <- as.dist(1 - as.matrix(sim))
n_tied <- sum(duplicated(as.numeric(d)))
msg("tied distances in the fixture similarity matrix:", n_tied,
    "of", length(d), "(expected: DreaMS scores are rounded to 3 decimals)")

msg("feature ids:", paste(mmo$feature_data$id, collapse = ", "))
msg("fixture sizes:")
print(file.size(list.files(FIXTURE, full.names = TRUE)) |>
        setNames(basename(list.files(FIXTURE))))
msg("total fixture bytes:", sum(file.size(list.files(FIXTURE, full.names = TRUE))))
msg("done")
