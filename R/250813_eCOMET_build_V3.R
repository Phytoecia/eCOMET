########################################################################################
# Define functions for data import and normalization
########################################################################################
# When this file is sourced directly (outside the package build), make sure the
# internal dependency checker exists so downstream calls to .require_pkg succeed.
if (!exists(".require_pkg", mode = "function")) {
  .require_pkg <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(pkg, " is required for this function. Install it with install.packages('", pkg, "').", call. = FALSE)
    }
  }
}

#' Import MZmine feature table and metadata to create a mmo (mass-spectrometry metabolomics object)
#' @description
#' Reads an MZmine exported feature table (typically the full feature table) and a sample metadata table
#' to initialize a mmo containing:
#' \itemize{
#'   \item \code{mmo$feature_data}: feature-by-sample abundance matrix (peak areas)
#'   \item \code{mmo$feature_info}: feature-level annotations (e.g., mz, rt, ranges, IDs)
#'   \item \code{mmo$metadata}: sample metadata with standardized \code{sample} and \code{group} columns
#' }
#' Sample columns in the MZmine table are matched to metadata using fuzzy string matching on area column names.
#'
#' @param mzmine_dir Path to the MZmine feature table CSV (should include feature-level columns and per-sample area columns)
#' @param metadata_dir Path to the metadata CSV file (must include sample_col and group_col)
#' @param group_col Column name in metadata used for grouping samples (e.g., treatment/species)
#' @param sample_col Column name in metadata used to identify and match samples to MZmine area columns
#' @param drop_missing_samples Logical. If FALSE (default), error when metadata samples are missing from
#' the MZmine table area columns. If TRUE, drop those samples from metadata (with a warning) and continue.
#' @param mz_col Optional m/z column name in the MZmine table (defaults to "mz" or "row m/z")
#' @param rt_col Optional RT column name in the MZmine table (defaults to "rt" or "row retention time")
#' @param max_distance Maximum edit distance used when fuzzy-matching metadata sample names to MZmine
#' area column names (default 5). Lower this for stricter matching.
#' @param feature_info_cols Character vector of feature-level columns to retain in \code{mmo$feature_info}.
#' Columns not present in the MZmine table are skipped with a warning.
#' @param pa_threshold Minimum peak area threshold for a feature to be considered present in a sample.
#' @return A mmo
#' @export
GetMZmineFeature <- function(mzmine_dir, metadata_dir, group_col, sample_col,
                             drop_missing_samples = FALSE,
                             mz_col = NULL, rt_col = NULL,
                             max_distance = 5,
                             pa_threshold = 1,
                             feature_info_cols = c(
                               "id",
                               "rt", "rt_range:min", "rt_range:max",
                               "mz", "mz_range:min", "mz_range:max",
                               "feature_group",
                               "ion_identities:iin_id",
                               "ion_identities:ion_identities"
                             )) {

  mmo <- list()

  # MZmine table: fread is 5-50x faster than read.csv on large files
  data <- data.table::fread(
    mzmine_dir,
    check.names  = FALSE,
    na.strings   = c("", "NA"),
    data.table   = FALSE
  )

  metadata <- data.table::fread(
    metadata_dir,
    check.names  = FALSE,
    na.strings   = c("", "NA"),
    data.table   = FALSE
  )
  # Check if group_col and sample_col are provided and exist in the metadata file
  if (missing(group_col) || !(group_col %in% names(metadata))) {
    stop("group_col must be provided and exist in the metadata file.", call. = FALSE)
  }
  if (missing(sample_col) || !(sample_col %in% names(metadata))) {
    stop("sample_col must be provided and exist in the metadata file.", call. = FALSE)
  }

  # --- detect mz / rt (or use overrides) ---
  if (is.null(mz_col)) {
    mz_col <- if ("mz" %in% names(data)) "mz"
    else if ("row m/z" %in% names(data)) "row m/z"
    else stop("No m/z column found (expected 'mz' or 'row m/z', or provide mz_col).", call. = FALSE)
  }
  if (is.null(rt_col)) {
    rt_col <- if ("rt" %in% names(data)) "rt"
    else if ("row retention time" %in% names(data)) "row retention time"
    else stop("No RT column found (expected 'rt' or 'row retention time', or provide rt_col).", call. = FALSE)
  }

  # --- generate stable feature key by mz_rt---
  data <- data |>
    dplyr::mutate(feature = paste(.data[[mz_col]], .data[[rt_col]], sep = "_"))

  # --- extract sample names from metadata ---
  samples <- trimws(as.character(metadata[[sample_col]]))
  samples_core <- sub("\\.(mzML|mzXML|raw)$", "", samples, ignore.case = TRUE)

  # Helper: fuzzy-match MZmine "area" columns to metadata sample names
  find_area_columns <- function(data, samples_core, samples_full, max_distance = 5) {
    colnames_all <- names(data)

    area_idx <- grepl("area", colnames_all, ignore.case = TRUE)
    if (!any(area_idx)) {
      stop("No columns in the MZmine table contain the word 'area' (case-insensitive).", call. = FALSE)
    }
    candidates <- colnames_all[area_idx]

    strip_wrappers <- function(x) {
      x <- sub("^datafile[:.]", "", x)                                        # datafile: or datafile.
      x <- sub("(?i)( Peak area|:area|\\.area)$", "", x, perl = TRUE)         # suffix variants
      x <- sub("(?i)\\.(mzml|mzxml|raw)$", "", x, perl = TRUE)                # trailing extensions
      x
    }
    cand_core <- strip_wrappers(candidates)

    D <- utils::adist(samples_core, cand_core, ignore.case = TRUE)
    best_idx <- apply(D, 1L, which.min)
    best_d   <- D[cbind(seq_along(samples_core), best_idx)]
    keep <- is.finite(best_d) & best_d <= max_distance

    list(
      area_cols = candidates[best_idx[keep]],
      keep_idx = keep,
      missing_full = samples_full[!keep]
    )
  }

  match_res <- find_area_columns(
    data = data,
    samples_core = samples_core,
    samples_full = samples,
    max_distance = max_distance
  )

  missing_samples_full <- match_res$missing_full
  keep_idx <- match_res$keep_idx
  area_columns <- match_res$area_cols

  if (length(area_columns) == 0) {
    stop(
      "No area columns matched any samples from metadata within max_distance.\n",
      "Check that you supplied an MZmine feature table containing per-sample area columns and that sample names correspond to metadata.",
      call. = FALSE
    )
  }

  if (any(duplicated(area_columns))) {
    dup <- unique(area_columns[duplicated(area_columns)])
    stop(
      "Multiple samples matched the same area column(s): ",
      paste(dup, collapse = ", "),
      ". Consider renaming columns or lowering max_distance.",
      call. = FALSE
    )
  }

  # --- handle missing samples (default ERROR with guidance; optional drop) ---
  if (length(missing_samples_full) > 0) {

    msg <- paste0(
      "Metadata contains ", length(missing_samples_full),
      " sample(s) that do not have matching 'area' columns in the MZmine table:\n",
      paste(missing_samples_full, collapse = ", "),
      "\n\nThis can occur if MZmine removed columns during processing (e.g., blank subtraction removing all features).\n\n",
      "To proceed anyway and drop these samples from metadata, re-run with `drop_missing_samples = TRUE`."
    )

    if (!isTRUE(drop_missing_samples)) {
      stop(msg, call. = FALSE)
    } else {
      warning(paste0(msg, "\n\nProceeding because `drop_missing_samples = TRUE`."), call. = FALSE)
    }
  }

  # Enforce alignment (either no missing, or user allowed dropping)
  metadata <- metadata[keep_idx, , drop = FALSE]
  samples <- samples[keep_idx]
  samples_core <- samples_core[keep_idx]

  # --- coerce area columns to numeric (handle rare comma-formatted exports) ---
  data[area_columns] <- lapply(data[area_columns], function(x) {
    if (is.character(x)) as.numeric(gsub(",", "", x)) else as.numeric(x)
  })

  # --- build feature abundance table (mmo$feature_data) ---
  if (!("id" %in% names(data))) {
    stop("No feature id column identified in the MZmine table (expected column 'id').", call. = FALSE)
  }

  feature_df <- data |>
    dplyr::select(.data$id, .data$feature, dplyr::all_of(area_columns)) |>
    dplyr::rename_with(~ samples_core, .cols = dplyr::all_of(area_columns))

  feature_df$id <- gsub(" ", "", as.character(feature_df$id))

  # --- Replace NA by 0 --- MC
  feature_df[is.na(feature_df)] <- 0

  # --- build feature info table (mmo$feature_info) ---
  # minimum required to be a useful feature_info
  required_min <- c("id", mz_col, rt_col)
  missing_min <- setdiff(required_min, names(data))
  if (length(missing_min) > 0) {
    stop(
      "The MZmine table is missing required feature-level columns: ",
      paste(missing_min, collapse = ", "),
      call. = FALSE
    )
  }

  # user-requested columns, but allow them to be absent (warn + skip)
  present_info <- intersect(feature_info_cols, names(data))
  absent_info  <- setdiff(feature_info_cols, names(data))
  if (length(absent_info) > 0) {
    warning(
      "Some requested feature_info columns were not present in the MZmine table and were skipped: ",
      paste(absent_info, collapse = ", "),
      call. = FALSE
    )
  }

  # Ensure key columns exist and are included
  present_info <- unique(c(present_info, "id", "feature"))

  feature_info <- data[, present_info, drop = FALSE]
  feature_info$id <- gsub(" ", "", as.character(feature_info$id))

  # Enforce column order: id, feature, then the rest
  col_order <- c("id", "feature", setdiff(names(feature_info), c("id", "feature")))
  feature_info <- feature_info[, col_order, drop = FALSE]
  rownames(feature_info) <- NULL

  # --- finalize metadata and output ---
  metadata$sample <- samples_core
  metadata$sample_full_exact <- samples
  metadata$group <- as.factor(metadata[[group_col]])

  mmo$feature_data <- feature_df
  mmo$feature_info <- feature_info
  mmo$pairwise <- data.frame(feature = mmo$feature_data$feature, id = mmo$feature_data$id) # make blank form for PairwiseComp
  mmo$metadata <- metadata

  class(mmo) <- "mmo"

  # Generate presence/absence table using set threshold
  mmo <- FeaturePresence(mmo, threshold = pa_threshold)

  message("Mass-spectrometry Metabolomics Object (mmo) created.")
  message(paste0("Feature number: ", nrow(mmo$feature_data)))
  message(paste0(nrow(mmo$metadata), " samples in ", length(unique(mmo$metadata$group)), " groups"))
  message(paste0("Presence/absence table generated with threshold: ", pa_threshold, " (PA)"))

  return(mmo)
}




#' Switch the group column in the \code{mmo}
#'
#' @description
#' This function switches the group column in the metadata of the \code{mmo} to a new specified column.
#' The new group column must exist in the metadata file.
#' @param mmo The \code{mmo}
#' @param new_group_col The name of the new group column in the metadata file
#' @return The \code{mmo} with the updated group column
#' @export
#' @examplesIf FALSE
#' mmo <- SwitchGroup(mmo, new_group_col = "genotype")
SwitchGroup <- function(mmo, new_group_col) {
  if (missing(new_group_col) || !(new_group_col %in% colnames(mmo$metadata))) {
    stop("new_group_col must be provided and must exist in the metadata file.")
  }
  mmo$metadata$group <- as.factor(mmo$metadata[[new_group_col]])
  print(paste0('Group column switched to ', new_group_col))
  print(paste0(length(unique(mmo$metadata$group)), ' groups in total'))
  print(paste0('The list of groups are: ', paste(levels(mmo$metadata$group), collapse = ', ')))
  return(mmo)
}

#' Adding annotation from SIRIUS to the \code{mmo}
#'
#' @description
#' This function reads SIRIUS structure identification and formula summary files,
#' and adds the annotations to the \code{mmo}.
#' @param mmo The \code{mmo}
#' @param canopus_structuredir Path to the SIRIUS structure_identification.tsv file
#' @param canopus_formuladir Path to the SIRIUS canopus_formula_summary.tsv file
#' @param filter_annot Logical. If TRUE, filter the annotations by probability threshold in CANOPUS.
#' @param filter_threshold Numeric between 0 and 1. The probability threshold for filtering annotations.
#' @return The \code{mmo} with SIRIUS annotations added
#' @export
#' @examplesIf FALSE
#' mmo <- AddSiriusAnnot(mmo,
#'  canopus_structuredir = "path/to/structure_identification.tsv",
#'  canopus_formuladir = "path/to/canopus_formula_summary.tsv"
#' )
#'
AddSiriusAnnot <- function(mmo, canopus_structuredir, canopus_formuladir, filter_annot = FALSE, filter_threshold = 0.5){
    structure_identifications <- readr::read_tsv(canopus_structuredir,
        show_col_types = FALSE)
    structure_identifications$mappingFeatureId <- gsub(" ", "",
        structure_identifications$mappingFeatureId)
    canopus_formula_summary <- readr::read_tsv(canopus_formuladir,
        show_col_types = FALSE)
    canopus_formula_summary$mappingFeatureId <- gsub(" ", "",
        canopus_formula_summary$mappingFeatureId)
    all_ids <- union(structure_identifications$mappingFeatureId,
        canopus_formula_summary$mappingFeatureId)
    duplicated_ids <- all_ids[duplicated(all_ids)]
    if (length(duplicated_ids) > 0) {
        warning(paste0("Duplicated IDs found in SIRIUS files: ",
            paste(duplicated_ids, collapse = ", ")))
    }
    sirius_df <- dplyr::select(mmo$feature_data, .data$id, .data$feature)
    structure_identifications$id <- structure_identifications$mappingFeatureId
    canopus_formula_summary$id <- canopus_formula_summary$mappingFeatureId

    shared_columns <- intersect(colnames(structure_identifications), colnames(canopus_formula_summary))
    shared_columns <- shared_columns[!shared_columns %in% c('id')]
    canopus_formula_summary <- canopus_formula_summary |> select(-all_of(shared_columns))

    sirius_df <- left_join(left_join(sirius_df, structure_identifications,
        by = c(id = "id"), multiple = "last"), canopus_formula_summary,
        by = 'id', multiple = "last")
    mmo$sirius_annot <- sirius_df
    print("SIRIUS ver6 is anticipated. Import may not be complete in lower versions.")
    print("SIRIUS annotation added to mmo$sirius_annot")
    if (filter_annot == TRUE){
      mmo <- filter_canopus_annotations(mmo, input = "sirius_annot", pathway_level = 'All', threshold = filter_threshold, na_prob = TRUE, overwrite = TRUE)
    }
    return(mmo)
}

#' Filter CANOPUS / SIRIUS annotations in an eCOMET \code{mmo} by probability threshold
#'
#' Applies a minimum-probability cutoff to selected CANOPUS (NPClassifier / ClassyFire)
#' annotation levels inside an eCOMET \code{mmo}. The function reads a chosen
#' annotation table from \code{mmo[[input]]} (default \code{"sirius_annot"}), flags
#' annotations below \code{threshold} by setting them to \code{NA}, and stores the result
#' as a new element on \code{mmo} named \code{"sirius_annot_filtered_<suffix>"}.
#'
#' Rows are never dropped. "Removing" an annotation means setting the annotation value
#' (and optionally its probability column) to \code{NA}.
#'
#' For each annotation column (e.g., \code{"NPC#pathway"}), the function looks for an
#' associated probability column using common SIRIUS export naming:
#' \itemize{
#'   \item \code{"<header> Probability"}
#'   \item \code{"<header> probability"}
#'   \item \code{"<header>Probability"}
#' }
#'
#' @param mmo An eCOMET \code{mmo} containing \code{mmo[[input]]} (a data.frame).
#' @param input Character. Name of the element on \code{mmo} to filter. Defaults to \code{"sirius_annot"}.
#'
#' @param pathway_level Character vector of one or more annotation header(s) to filter.
#'   Valid options include:
#'   \itemize{
#'     \item \code{"All"}
#'     \item \code{"All_NPC"}
#'     \item \code{"All_ClassyFire"}
#'     \item \code{"NPC#pathway"}
#'     \item \code{"NPC#superclass"}
#'     \item \code{"NPC#class"}
#'     \item \code{"ClassyFire#superclass"}
#'     \item \code{"ClassyFire#class"}
#'     \item \code{"ClassyFire#subclass"}
#'     \item \code{"ClassyFire#level 5"}
#'     \item \code{"ClassyFire#most specific class"}
#'   }
#'   Special values:
#'   \itemize{
#'     \item \code{"All"} expands to all NPC + ClassyFire levels listed above.
#'     \item \code{"All_NPC"} expands to NPC levels only.
#'     \item \code{"All_ClassyFire"} expands to ClassyFire levels only.
#'   }
#'
#' @param threshold Decimal between 0 and 1. Annotations with probability < threshold are flagged to \code{NA}.
#' @param na_prob Logical. If \code{TRUE} (default), also set the corresponding probability value to \code{NA}.
#'
#' @param suffix Optional character string appended to the created element name:
#'   \code{"sirius_annot_filtered_<suffix>"}. If \code{NULL} (default), a suffix is auto-generated.
#'
#' @param overwrite Logical. If \code{FALSE} (default) and the target element already exists on \code{mmo},
#'   the function errors to avoid accidental overwrites.
#'
#' @param verbose Logical. If \code{TRUE}, prints a concise summary including counts of
#'   non-missing annotations before and after filtering.
#'
#' @return The updated \code{mmo}, with a new element \code{mmo[[paste0("sirius_annot_filtered_", suffix)]]}.
#'
#' @examples
#' \dontrun{
#' mmo <- filter_canopus_annotations(mmo, input = "sirius_annot",
#'                                  pathway_level = "NPC#pathway", threshold = 0.9,
#'                                  suffix = "NPC_pathway_0.9", verbose = TRUE)
#'
#' mmo <- filter_canopus_annotations(mmo, input = "sirius_annot_filtered_NPC_pathway_0.9",
#'                                  pathway_level = "All_NPC", threshold = 0.95,
#'                                  suffix = "All_NPC_0.95", verbose = TRUE)
#' }
#'
#' @export
filter_canopus_annotations <- function(
    mmo,
    input = "sirius_annot",
    pathway_level = "NPC#pathway",
    threshold = 0.9,
    na_prob = TRUE,
    suffix = NULL,
    overwrite = FALSE,
    verbose = TRUE
) {
  # Basic checks
  if (is.null(mmo) || !is.list(mmo)) {
    stop("`mmo` must be a list-like ecomet object.", call. = FALSE)
  }
  if (!is.character(input) || length(input) != 1L || is.na(input) || !nzchar(input)) {
    stop("`input` must be a single, non-empty character string.", call. = FALSE)
  }
  if (!(input %in% names(mmo)) || !is.data.frame(mmo[[input]])) {
    stop("`mmo[['", input, "']]` must exist and be a data.frame.", call. = FALSE)
  }
  if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold) ||
      threshold < 0 || threshold > 1) {
    stop("`threshold` must be a single numeric value in [0, 1].", call. = FALSE)
  }

  df <- mmo[[input]]

  # Allowed levels and expansions
  npc_levels <- c("NPC#pathway", "NPC#superclass", "NPC#class")
  cf_levels  <- c("ClassyFire#superclass", "ClassyFire#class", "ClassyFire#subclass",
                  "ClassyFire#level 5", "ClassyFire#most specific class")
  all_levels <- c(npc_levels, cf_levels)
  allowed <- c("All", "All_NPC", "All_ClassyFire", all_levels)

  if (is.list(pathway_level)) pathway_level <- unlist(pathway_level, use.names = FALSE)
  pathway_level <- as.character(pathway_level)

  # Expand selectors
  if (any(pathway_level %in% "All")) {
    pathway_level <- unique(c(setdiff(pathway_level, "All"), all_levels))
  }
  if (any(pathway_level %in% "All_NPC")) {
    pathway_level <- unique(c(setdiff(pathway_level, "All_NPC"), npc_levels))
  }
  if (any(pathway_level %in% "All_ClassyFire")) {
    pathway_level <- unique(c(setdiff(pathway_level, "All_ClassyFire"), cf_levels))
  }

  bad <- setdiff(pathway_level, allowed)
  if (length(bad) > 0) {
    stop("Unknown `pathway_level` value(s): ", paste0(bad, collapse = ", "), call. = FALSE)
  }

  # Helper: find probability column variants for a given header
  find_prob_col <- function(df, header) {
    candidates <- c(
      paste0(header, " Probability"),
      paste0(header, " probability"),
      paste0(header, "Probability")
    )
    hit <- candidates[candidates %in% names(df)]
    if (length(hit) == 0) return(NA_character_)
    hit[[1]]
  }

  # Only operate on headers that exist as annotation columns
  present_headers <- pathway_level[pathway_level %in% names(df)]
  if (length(present_headers) == 0) {
    stop(
      "None of the requested `pathway_level` columns were found in `mmo[['", input, "']]`.\n",
      "Requested: ", paste(pathway_level, collapse = ", "),
      call. = FALSE
    )
  }

  # Pair each header with a probability column; skip levels without prob col
  prob_cols <- vapply(present_headers, function(h) find_prob_col(df, h), character(1))
  no_prob <- present_headers[is.na(prob_cols)]
  if (length(no_prob) > 0) {
    warning(
      "Probability column not found for: ",
      paste0(no_prob, collapse = ", "),
      "\nExpected e.g. '<header> Probability'. These levels will be skipped.",
      call. = FALSE
    )
    keep_idx <- !is.na(prob_cols)
    present_headers <- present_headers[keep_idx]
    prob_cols <- prob_cols[keep_idx]
  }
  if (length(present_headers) == 0) {
    stop("No usable (header + probability) pairs were found for filtering.", call. = FALSE)
  }

  # Auto suffix if not provided (sanitize for safe-ish $ access)
  if (is.null(suffix) || is.na(suffix) || !nzchar(suffix)) {
    base <- pathway_level[[1]]
    base <- gsub("#", "_", base, fixed = TRUE)
    base <- gsub("[^A-Za-z0-9_\\.]+", "_", base)
    thr  <- formatC(threshold, format = "f", digits = 2)
    suffix <- paste0(base, "_", thr)
  }

  target_name <- paste0("sirius_annot_filtered_", suffix)
  if (!overwrite && (target_name %in% names(mmo))) {
    stop(
      "An element named '", target_name, "' already exists on `mmo`.\n",
      "Set `overwrite = TRUE` or choose a different `suffix`.",
      call. = FALSE
    )
  }

  # Count non-missing annotations BEFORE (across all selected levels)
  count_nonmissing <- function(vec) {
    v <- vec
    if (!is.character(v)) v <- as.character(v)
    sum(!is.na(v) & nzchar(v))
  }
  before_count <- sum(vapply(present_headers, function(h) count_nonmissing(df[[h]]), integer(1)))

  # Filter: set annotations (and optionally probs) to NA where prob < threshold OR missing prob
  out <- df

  removed_total <- 0L
  for (j in seq_along(present_headers)) {
    h <- present_headers[[j]]
    pcol <- prob_cols[[j]]

    p <- suppressWarnings(as.numeric(out[[pcol]]))
    ann <- out[[h]]
    ann_chr <- if (is.character(ann)) ann else as.character(ann)

    failing <- is.na(ann_chr) | !nzchar(ann_chr) | is.na(p) | (p < threshold)

    was_present <- !is.na(ann_chr) & nzchar(ann_chr)
    removed_here <- sum(was_present & failing)
    removed_total <- removed_total + removed_here

    out[[h]][failing] <- NA
    if (isTRUE(na_prob)) out[[pcol]][failing] <- NA_real_
  }

  # Count non-missing annotations AFTER (across all selected levels)
  after_count <- sum(vapply(present_headers, function(h) count_nonmissing(out[[h]]), integer(1)))

  mmo[[target_name]] <- out

  if (verbose) {
    message(
      "filter_canopus_annotations(): stored as mmo[['", target_name, "']]\n",
      "Input: mmo[['", input, "']]\n",
      "Levels filtered: ", paste(present_headers, collapse = ", "), "\n",
      "Threshold (kept >=): ", threshold, "\n",
      "Non-missing annotations (before -> after): ", before_count, " -> ", after_count, "\n",
      "Removed annotations below threshold: ", removed_total
    )
  }

  mmo
}



#' Filter SIRIUS structure (CSI:FingerID) annotations by COSMIC confidence score
#'
#' Applies a minimum COSMIC confidence cutoff to SIRIUS structure predictions inside an
#' ecomet \code{mmo}. The function reads a chosen annotation table from
#' \code{mmo[[input]]} (default \code{"sirius_annot"}), flags structure annotations below
#' \code{threshold} by setting selected structure-identification fields to \code{NA},
#' and stores the result as a new element on \code{mmo} named
#' \code{"sirius_annot_filtered_<suffix>"}.
#'
#' Rows are never dropped. "Removing" a structure means setting selected fields
#' (e.g., SMILES/InChI/name/InChIkey2D) to \code{NA}.
#'
#' @param mmo An eCOMET \code{mmo} containing \code{mmo[[input]]} (a data.frame).
#' @param input Character. Name of the element on \code{mmo} to filter. Defaults to \code{"sirius_annot"}.
#'
#' @param cosmic_mode Which COSMIC column to use. One of \code{"exact"} or \code{"approx"}.
#' @param threshold Numeric. Keep structures with COSMIC >= threshold.
#'
#' @param fields Character vector of columns to NA-out when COSMIC < threshold.
#'   If \code{"auto"} (default), uses a reasonable default set if present in the table.
#'
#' @param na_cosmic Logical. If \code{TRUE} (default), also set the COSMIC value to \code{NA}
#'   when the structure is filtered out.
#'
#' @param suffix Optional character string appended to the created element name.
#' @param overwrite Logical. If \code{FALSE} (default) and the target element already exists, error.
#' @param verbose Logical. If \code{TRUE}, prints a concise summary.
#'
#' @return The updated \code{mmo}, with a new element \code{mmo[[paste0("sirius_annot_filtered_", suffix)]]}.
#'
#' @examples
#' \dontrun{
#' mmo <- filter_cosmic_structure(mmo, input = "sirius_annot",
#'                               cosmic_mode = "approx", threshold = 0.5,
#'                               suffix = "COSMIC_exact_0.5", verbose = TRUE)
#' }
#'
#' @export
filter_cosmic_structure <- function(
    mmo,
    input = "sirius_annot",
    cosmic_mode = c("exact", "approx"),
    threshold = 0.5,
    fields = "auto",
    na_cosmic = TRUE,
    suffix = NULL,
    overwrite = FALSE,
    verbose = TRUE
) {
  # Basic checks
  if (is.null(mmo) || !is.list(mmo)) {
    stop("`mmo` must be a list-like ecomet object.", call. = FALSE)
  }
  if (!is.character(input) || length(input) != 1L || is.na(input) || !nzchar(input)) {
    stop("`input` must be a single, non-empty character string.", call. = FALSE)
  }
  if (!(input %in% names(mmo)) || !is.data.frame(mmo[[input]])) {
    stop("`mmo[['", input, "']]` must exist and be a data.frame.", call. = FALSE)
  }
  cosmic_mode <- match.arg(cosmic_mode)

  if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold)) {
    stop("`threshold` must be a single non-missing numeric value.", call. = FALSE)
  }

  df <- mmo[[input]]

  cosmic_col <- switch(
    cosmic_mode,
    exact  = "ConfidenceScoreExact",
    approx = "ConfidenceScoreApproximate"
  )

  if (!(cosmic_col %in% names(df))) {
    stop(
      "COSMIC column '", cosmic_col, "' not found in `mmo[['", input, "']]`.",
      call. = FALSE
    )
  }

  # Determine fields to blank when failing
  if (length(fields) == 1L && identical(fields, "auto")) {
    candidate_fields <- c(
      # core structure / identifiers often present in SIRIUS exports
      "smiles", "InChI", "InChIkey2D", "name",
      # database linkage / provenance
      "links", "pubchemids", "dbflags",
      # ranking/scoring columns (structure side)
      "structurePerIdRank", "CSI:FingerIDScore",
      # optional: formula/overall scores that are often carried along
      "SiriusScore", "SiriusScoreNormalized",
      # COSMIC columns (optionally NA'd below)
      "ConfidenceScoreExact", "ConfidenceScoreApproximate"
    )
    fields_to_na <- intersect(candidate_fields, names(df))
  } else {
    fields <- as.character(fields)
    fields_to_na <- intersect(fields, names(df))
    if (length(fields_to_na) == 0) {
      stop("None of the requested `fields` were found in `mmo[['", input, "']]`.", call. = FALSE)
    }
  }

  # Helper: count "non-missing structure annotation" for reporting
  count_nonmissing <- function(vec) {
    v <- vec
    if (!is.character(v)) v <- as.character(v)
    sum(!is.na(v) & nzchar(v))
  }

  # Choose an indicator column for "structure exists" (best-effort)
  indicator_col <- if ("smiles" %in% names(df)) {
    "smiles"
  } else if ("InChI" %in% names(df)) {
    "InChI"
  } else if ("name" %in% names(df)) {
    "name"
  } else {
    cosmic_col
  }

  before_count <- count_nonmissing(df[[indicator_col]])

  out <- df

  # Parse COSMIC, handling character, NA, and +/-Inf
  cosmic <- suppressWarnings(as.numeric(out[[cosmic_col]]))

  ind <- out[[indicator_col]]
  ind_chr <- if (is.character(ind)) ind else as.character(ind)
  has_struct <- !is.na(ind_chr) & nzchar(ind_chr)

  # Fail if structure present but COSMIC missing/non-finite/below threshold
  failing <- has_struct & (is.na(cosmic) | !is.finite(cosmic) | (cosmic < threshold))

  removed_total <- sum(failing, na.rm = TRUE)

  if (removed_total > 0) {
    for (col in fields_to_na) {
      out[[col]][failing] <- NA
    }
    if (isTRUE(na_cosmic)) {
      out[[cosmic_col]][failing] <- NA_real_
    }
  }

  after_count <- count_nonmissing(out[[indicator_col]])

  # Auto suffix if not provided
  if (is.null(suffix) || is.na(suffix) || !nzchar(suffix)) {
    thr <- formatC(threshold, format = "f", digits = 2)
    suffix <- paste0("COSMIC_", cosmic_mode, "_", thr)
  }
  suffix <- gsub("[^A-Za-z0-9_\\.]+", "_", suffix)

  target_name <- paste0("sirius_annot_filtered_", suffix)
  if (!overwrite && (target_name %in% names(mmo))) {
    stop(
      "An element named '", target_name, "' already exists on `mmo`.\n",
      "Set `overwrite = TRUE` or choose a different `suffix`.",
      call. = FALSE
    )
  }

  mmo[[target_name]] <- out

  if (verbose) {
    message(
      "filter_cosmic_structure(): stored as mmo[['", target_name, "']]\n",
      "Input: mmo[['", input, "']]\n",
      "COSMIC mode: ", cosmic_mode, " (column: ", cosmic_col, ")\n",
      "Threshold (kept >=): ", threshold, "\n",
      "Structure indicator column: ", indicator_col, "\n",
      "Non-missing structures (before -> after): ", before_count, " -> ", after_count, "\n",
      "Removed structures below threshold: ", removed_total, "\n",
      "Fields NA'd: ", paste(fields_to_na, collapse = ", ")
    )
  }

  mmo
}



#' Add custom annotations to a \code{mmo}
#'
#' @description
#' Match features to a custom DB by m/z (ppm) and RT (minutes) tolerances and
#' attach a list-column of candidate compound IDs per feature.
#'
#' @param mmo A \code{mmo} created by `GetMZmineFeature()`.
#' @param DB_file CSV path with at least columns `compound`, `mz`, `rt`.
#' @param mztol m/z tolerance in ppm (default 5).
#' @param rttol RT tolerance in minutes (default 0.5).
#' @return The same `mmo` with `mmo$custom_annot` (id, feature, custom_annot).
#' @export
#' @examplesIf FALSE
#' mmo <- AddCustomAnnot(mmo, DB_file = "path/to/custom_db.csv", mztol = 5, rttol = 0.5)
AddCustomAnnot <- function(mmo, DB_file, mztol = 5, rttol = 0.5) {
  DB <- read.csv(DB_file, stringsAsFactors = FALSE)

  DB <- dplyr::mutate(
    DB,
    mz = as.numeric(.data$mz),
    rt = as.numeric(.data$rt)
  )

  # Parse "feature" (formatted like "mz_rt") into numeric mz/rt columns
  feature_annot <- mmo$feature_data |>
    tidyr::separate(
      col = "feature",
      into = c("mz", "rt"),
      sep = "_",
      remove = FALSE,
      convert = TRUE
    )

  # For each feature, collect matching DB compounds within tolerances
  annotated_features <- dplyr::mutate(
    feature_annot,
    custom_annot = purrr::map2(.data$mz, .data$rt, function(mzi, rti) {
      dplyr::filter(
        DB,
        abs(mzi - .data$mz) / mzi * 1e6 <= mztol,
        abs(rti - .data$rt) <= rttol
      ) |>
        dplyr::pull(.data$compound)
    })
  )

  mmo$custom_annot <- annotated_features |> dplyr::select(.data$id, .data$feature, .data$custom_annot)

  message("Custom annotation added to mmo$custom_annot using ", DB_file)
  mmo
}




#' #' Replace zero and NA values in the \code{mmo}
#'
#' This function replaces zero values in the feature table of a \code{mmo}.
#' Imputed data are stored in mmo$imputed_feature_data, to be used for
#' downsteam analyses including PairwiseComp().
#' Note that imputation affects normalizations (Log-transformation, etc.),
#' as well as chemical diversity calculations that uses presence/absence.
#'
#' @param mmo A \code{mmo} containing `feature_data`
#' @param method Replacement method:
#'   - "one": replace zeros and NA values with 1
#'   - "half_min": replace zeros and NA values with half of the smallest
#'     non-zero value in each feature (row)
#'
#' @return The updated \code{mmo}
#' @export
ReplaceZero <- function(mmo, method = c("one", "half_min")) {
  df <- mmo$feature_data

  # assume first two columns are identifiers
  id_feat <- df[, 1:2, drop = FALSE]
  num_df  <- df[, -(1:2), drop = FALSE]

  # coerce numeric values only
  num_mat <- as.matrix(
    sapply(num_df, function(x) as.numeric(as.character(x)))
  )

  if (method == "one") {
    new_mat <- num_mat
    new_mat[is.na(new_mat) | new_mat == 0] <- 1L
  } else {
    # half_min: replace 0/NA with half the per-feature minimum positive value
    pos_mat <- num_mat
    pos_mat[is.na(pos_mat) | pos_mat <= 0] <- NA_real_
    row_mins <- apply(pos_mat, 1, min, na.rm = TRUE)
    if (any(is.infinite(row_mins))) {
      bad <- which(is.infinite(row_mins))
      stop("half_min replacement failed: feature(s) at row(s) ",
           paste(bad, collapse = ", "), " have no non-zero values.", call. = FALSE)
    }
    half_mins  <- row_mins / 2
    replace_at <- is.na(num_mat) | num_mat == 0
    new_mat    <- num_mat
    # broadcast per-row half_min to full matrix, then apply mask
    new_mat[replace_at] <- matrix(half_mins, nrow = nrow(num_mat),
                                  ncol = ncol(num_mat))[replace_at]
  }

  num_df_out <- as.data.frame(new_mat, stringsAsFactors = FALSE)
  names(num_df_out) <- names(num_df)

  mmo$imputed_feature_data <- cbind(id_feat, num_df_out)
  message(sprintf("0 values were replaced with %s, then stored as mmo$imputed_feature_data", method))
  mmo
}

#' Convert feature abundances to presence / absence
#'
#' This function converts the feature abundance matrix in a \code{mmo}
#' into a binary presence/absence matrix and stores it as a new component
#' of the \code{mmo} (\code{mmo$feature_presence}).
#'
#' A feature is considered present (1) if its abundance is greater than
#' a specified threshold, and absent (0) otherwise.
#'
#' This function does NOT overwrite \code{mmo$feature_data}.
#'
#' @param mmo The \code{mmo}
#' @param threshold Numeric threshold for presence (default = 1).
#'   Values > threshold are set to 1, values <= threshold or NA are set to 0.
#'
#' @return The \code{mmo} with a new presence/absence table
#'   stored in \code{mmo$feature_presence}
#'
#' @export
#'
#' @examplesIf FALSE
#' mmo <- FeaturePresence(mmo, threshold = 1)
FeaturePresence <- function(mmo, threshold = 1) {

  df <- mmo$feature_data

  # assume first two columns are id and feature
  id_feat <- df[, 1:2]
  num_df  <- df[, -(1:2), drop = FALSE]

  # force numeric matrix
  num_mat <- as.matrix(sapply(num_df, as.numeric))

  # convert to presence / absence
  pa_mat <- ifelse(num_mat > threshold, 1, 0)

  # rebuild data.frame
  pa_df <- as.data.frame(pa_mat, stringsAsFactors = FALSE)
  names(pa_df) <- names(num_df)

  df_out <- cbind(id_feat, pa_df)

  # store as new slot
  mmo$feature_presence <- df_out

  message(
    sprintf(
      "Feature presence/absence matrix created using threshold > %s",
      threshold
    )
  )

  return(mmo)
}



#' Use sample mass in the metadata file to normalize the peak area
#'
#' This function normalizes the peak area in the feature data of the \code{mmo} by the mass of each sample, provided in the metadata.
#' The feature data is replaced by (original value * mean mass) / sample mass.
#'
#' @param mmo The \code{mmo}
#' @return The \code{mmo} with normalized feature data (mmo$feature_data)
#' @export
#' @examplesIf FALSE
#' mmo <- MassNormalization(mmo)
MassNormalization <- function(mmo){
  metadata  <- mmo$metadata
  mean_mass <- mean(metadata$mass)

  id_cols   <- mmo$feature_data[, 1:2, drop = FALSE]
  mat       <- as.matrix(mmo$feature_data[, -(1:2), drop = FALSE])

  # Match sample order in the matrix to metadata$mass in one vectorised step
  mass_vec  <- metadata$mass[match(colnames(mat), metadata$sample)]
  mat       <- sweep(mat, 2, mean_mass / mass_vec, `*`)

  mmo$feature_data <- cbind(id_cols, as.data.frame(mat, check.names = FALSE))
  message("Peak area normalised by sample mass.")
  return(mmo)
}

#' Log-normalize the peak area in the \code{mmo}
#'
#' This function applies log2 transformation to the peak area in the feature data of the \code{mmo}.
#' Run ReplaceZero() before this function to avoid -Inf values.
#'
#' @param mmo The \code{mmo}
#' @param imputed_data Whether to use imputed feature data (default = FALSE)
#' @return The \code{mmo} with log-normalized feature data (mmo$log)
#' @export
#' @examplesIf FALSE
#' mmo <- LogNormalization(mmo)
LogNormalization <- function(mmo, imputed_data = FALSE){
  if (imputed_data == FALSE){
    feature_data_only <- mmo$feature_data[,-(1:2)]
    feature_ids <- mmo$feature_data[, 1:2]
  }else{
    feature_data_only <- mmo$imputed_feature_data[,-(1:2)]
    feature_ids <- mmo$imputed_feature_data[, 1:2]
  }
  log_data <- log2(feature_data_only+1)
  log_df <- cbind(feature_ids, log_data)
  mmo$log <- log_df
  print('Log-normalized values were added to mmo$log')
  return(mmo)
}

#' Mean-center the peak area in the \code{mmo}
#'
#' This function applies mean-centering to the peak area in the feature data
#' of the \code{mmo}. Mean-centering is performed per feature (row) across samples.
#' Features with zero variance are returned as all zeros and are reported in a warning.
#'
#' @param mmo The \code{mmo}
#' @param imputed_data Whether to use imputed feature data (default = FALSE)
#' @return The \code{mmo} with mean-centered feature data stored in `mmo$meancentered`
#' @export
MeancenterNormalization <- function(mmo, imputed_data = FALSE){
  if (imputed_data == FALSE){
    feature_data_only <- mmo$feature_data[, -(1:2)]
  }else{
    feature_data_only <- mmo$imputed_feature_data[,-(1:2)]
  }
  feature_ids <- mmo$feature_data[, 1:2]

  # detect zero-variance features
  row_sd <- apply(feature_data_only, 1, sd, na.rm = TRUE)
  zero_sd_idx <- which(is.na(row_sd) | row_sd == 0)

  if (length(zero_sd_idx) > 0) {
    warning(
      "MeancenterNormalization: the following features have zero variance and were mean-centered to all zeros:\n",
      paste(
        apply(feature_ids[zero_sd_idx, , drop = FALSE], 1, paste, collapse = " | "),
        collapse = "\n"
      ),
      call. = FALSE
    )
  }

  mat <- as.matrix(feature_data_only)
  mean_centered_data <- sweep(mat, 1, rowMeans(mat, na.rm = TRUE), `-`)

  mean_centered_df <- cbind(feature_ids, mean_centered_data)
  mmo$meancentered <- mean_centered_df

  return(mmo)
}

#' Z-normalize the peak area in the \code{mmo}
#'
#' This function applies Z-score normalization to the peak area in the feature data
#' of the \code{mmo}. Z-scores are calculated per feature (row) across samples.
#' Features with zero variance cannot be Z-normalized and are returned as NA.
#'
#' @param mmo The \code{mmo}
#' @param imputed_data Whether to use imputed feature data (default = FALSE)
#' @return The \code{mmo} with Z-normalized feature data stored in `mmo$zscore`
#' @export
ZNormalization <- function(mmo, imputed_data = FALSE) {
  if (imputed_data == FALSE){
    feature_data_only <- mmo$feature_data[, -(1:2)]
  }else{
    feature_data_only <- mmo$imputed_feature_data[,-(1:2)]
  }
  feature_ids <- mmo$feature_data[, 1:2]

  # compute per-feature SD
  row_sd <- apply(feature_data_only, 1, sd, na.rm = TRUE)

  # identify zero-variance features
  zero_sd_idx <- which(is.na(row_sd) | row_sd == 0)

  if (length(zero_sd_idx) > 0) {
    warning(
      "ZNormalization: the following features have zero variance and were set to NA:\n",
      paste(
        apply(feature_ids[zero_sd_idx, , drop = FALSE], 1, paste, collapse = " | "),
        collapse = "\n"
      ),
      call. = FALSE
    )
  }

  mat        <- as.matrix(feature_data_only)
  row_means  <- rowMeans(mat, na.rm = TRUE)
  mat_c      <- sweep(mat, 1, row_means, `-`)
  zscore_mat <- sweep(mat_c, 1, row_sd, `/`)
  zscore_mat[is.na(row_sd) | row_sd == 0, ] <- NA_real_

  zscore_df <- cbind(feature_ids, zscore_mat)
  colnames(zscore_df) <- colnames(mmo$feature_data)
  mmo$zscore <- zscore_df

  return(mmo)
}


# DEPRECATED: superseded by AddChemSim().
# AddChemDist() stores 1-similarity (dissimilarity) in .dissim slots.
# AddChemSim() stores raw similarity in .sim slots and supports
# sparse matrices for large datasets. Retained for backwards
# compatibility until downstream functions are updated to .sim slots.
#' Add chemical distance matrices to the \code{mmo}
#'
#' This function reads cosine, DREAMS, and MS2DeepScore molecular networking outputs from MZmine,
#' then transform the similarity to distance and adds the dissimilarity matrices to the \code{mmo}.
#'
#' @param mmo The \code{mmo}
#' @param cos_dir Path to the cosine similarity CSV file from MZMine (molecular networking)
#' @param dreams_dir Path to the DREAMS similarity CSV file from MZMine (molecular networking)
#' @param m2ds_dir Path to the MS2DeepScore similarity CSV file from MZMine (molecular networking)
#' @return The \code{mmo} with dissimilarity matrices added (mmo$cos.dissim, mmo$dreams.dissim, mmo$m2ds.dissim)
#' @export
#' @examplesIf FALSE
#' mmo <- AddChemDist(mmo,
#'  cos_dir = "path/to/cosine_similarity.csv",
#'  dreams_dir = "path/to/dreams_similarity.csv",
#'  m2ds_dir = "path/to/ms2deepscore_similarity.csv"
#' )
AddChemDist <- function(mmo, cos_dir = NULL, dreams_dir = NULL, m2ds_dir = NULL) {
  .require_pkg("data.table")
  add_dissim_matrix <- function(mmo, sim_dir, slot_name) {
    # Fast read
    sim <- data.table::fread(sim_dir, col.names = c("cluster1", "cluster2", "metric", "similarity", "etc"))
    sim$dissimilarity <- 1 - sim$similarity

    # Ensure cluster IDs are clean characters
    sim$cluster1 <- trimws(as.character(sim$cluster1))
    sim$cluster2 <- trimws(as.character(sim$cluster2))

    # Build mapping
    clusters <- unique(c(sim$cluster1, sim$cluster2))
    cluster_index <- stats::setNames(seq_along(clusters), clusters)

    # Map clusters to integer indices
    i <- cluster_index[sim$cluster1]
    j <- cluster_index[sim$cluster2]

    # Check for problems
    if (anyNA(i) || anyNA(j)) {
      bad <- unique(c(sim$cluster1[is.na(i)], sim$cluster2[is.na(j)]))
      stop("Cluster IDs not found in mapping: ", paste(bad, collapse = ", "))
    }

    # Preallocate full dense matrix
    n <- length(clusters)
    dissim_mat <- matrix(1, nrow = n, ncol = n)
    diag(dissim_mat) <- 0

    # Fill dissimilarities
    dissim_mat[cbind(i, j)] <- sim$dissimilarity
    dissim_mat[cbind(j, i)] <- sim$dissimilarity

    # Add names
    dimnames(dissim_mat) <- list(clusters, clusters)

    # Reduce memory footprint
    mode(dissim_mat) <- "single"

    mmo[[slot_name]] <- dissim_mat
    message(slot_name, " added to mmo")
    return(mmo)
  }

  if (!is.null(cos_dir))    mmo <- add_dissim_matrix(mmo, cos_dir,    "cos.dissim")
  if (!is.null(dreams_dir)) mmo <- add_dissim_matrix(mmo, dreams_dir, "dreams.dissim")
  if (!is.null(m2ds_dir))   mmo <- add_dissim_matrix(mmo, m2ds_dir,   "m2ds.dissim")

  if (is.null(cos_dir) && is.null(dreams_dir) && is.null(m2ds_dir)) {
    stop("Please provide at least one valid directory.")
  }

  return(mmo)
}


#' Add chemical similarity matrices to the \code{mmo}
#'
#' This function reads cosine, DREAMS, and MS2DeepScore molecular networking outputs from MZmine
#' and stores the raw pairwise \strong{similarity} matrices (values 0--1, higher = more similar)
#' in the \code{mmo}. For datasets with more than 10,000 features a sparse \code{Matrix} is used
#' to avoid allocating a full dense matrix; unobserved feature pairs carry implicit similarity 0
#' (equivalent to dissimilarity 1, maximally different).
#'
#' @param mmo The \code{mmo}
#' @param cos_dir Path to the cosine similarity CSV file from MZMine (molecular networking)
#' @param dreams_dir Path to the DREAMS similarity CSV file from MZMine (molecular networking)
#' @param m2ds_dir Path to the MS2DeepScore similarity CSV file from MZMine (molecular networking)
#' @return The \code{mmo} with similarity matrices stored in \code{mmo$cos.sim},
#'   \code{mmo$dreams.sim}, and/or \code{mmo$m2ds.sim}. Large datasets (n > 10,000)
#'   use sparse \code{dgCMatrix} format; smaller datasets use a dense numeric matrix.
#' @note Tree-based downstream functions (\code{GetFaithPD}, \code{FeatureDendrogram},
#'   \code{GetBetaDiversity} with \code{method = "Gen.Uni"}) require a full n x n distance
#'   matrix and will fail at large scale regardless of sparse storage. Use \code{filter_mmo()}
#'   to reduce feature count before calling those functions.
#' @export
#' @examplesIf FALSE
#' mmo <- AddChemSim(mmo,
#'  cos_dir = "path/to/cosine_similarity.csv",
#'  dreams_dir = "path/to/dreams_similarity.csv",
#'  m2ds_dir = "path/to/ms2deepscore_similarity.csv"
#' )
AddChemSim <- function(mmo, cos_dir = NULL, dreams_dir = NULL, m2ds_dir = NULL) {
  .require_pkg("data.table")
  add_sim_matrix <- function(mmo, sim_dir, slot_name) {
    # Fast read
    sim <- data.table::fread(sim_dir, col.names = c("cluster1", "cluster2", "metric", "similarity", "etc"))

    # Ensure cluster IDs are clean characters
    sim$cluster1 <- trimws(as.character(sim$cluster1))
    sim$cluster2 <- trimws(as.character(sim$cluster2))

    # Build mapping
    clusters <- unique(c(sim$cluster1, sim$cluster2))
    cluster_index <- stats::setNames(seq_along(clusters), clusters)

    # Map clusters to integer indices
    i <- cluster_index[sim$cluster1]
    j <- cluster_index[sim$cluster2]

    # Check for problems
    if (anyNA(i) || anyNA(j)) {
      bad <- unique(c(sim$cluster1[is.na(i)], sim$cluster2[is.na(j)]))
      stop("Cluster IDs not found in mapping: ", paste(bad, collapse = ", "))
    }

    n <- length(clusters)

    if (n > 10000L) {
      message(
        "Large dataset (n = ", n, " features): similarity matrix stored in sparse format to save memory.\n",
        slot_name, " added to mmo\n\n",
        "Tree-based functions at this scale:\n",
        "  - GetFaithPD and GetBetaDiversity(method='Gen.Uni') will automatically use\n",
        "    minimum spanning tree (single-linkage) instead of average-linkage clustering.\n",
        "    Results are valid but may differ from smaller datasets.\n",
        "  - FeatureDendrogram requires a dense matrix and will not run at this scale.\n",
        "    Use filter_mmo() to subset features first."
      )

      # Include diagonal (self-similarity = 1) in the sparse triplets
      diag_idx <- seq_len(n)
      i_all   <- c(i, j, diag_idx)
      j_all   <- c(j, i, diag_idx)
      x_all   <- c(sim$similarity, sim$similarity, rep(1, n))

      # Build sparse similarity matrix; implicit 0 = similarity 0 = dissimilarity 1
      sim_mat <- Matrix::sparseMatrix(
        i = i_all, j = j_all, x = x_all,
        dims = c(n, n),
        dimnames = list(clusters, clusters)
      )
    } else {
      # Preallocate full dense matrix; 0 = no observed similarity, diagonal = self-similarity 1
      sim_mat <- matrix(0, nrow = n, ncol = n)
      diag(sim_mat) <- 1

      # Fill similarities
      sim_mat[cbind(i, j)] <- sim$similarity
      sim_mat[cbind(j, i)] <- sim$similarity

      # Add names
      dimnames(sim_mat) <- list(clusters, clusters)

      # Reduce memory footprint
      mode(sim_mat) <- "single"
    }

    mmo[[slot_name]] <- sim_mat
    message(slot_name, " added to mmo")
    return(mmo)
  }

  if (!is.null(cos_dir))    mmo <- add_sim_matrix(mmo, cos_dir,    "cos.sim")
  if (!is.null(dreams_dir)) mmo <- add_sim_matrix(mmo, dreams_dir, "dreams.sim")
  if (!is.null(m2ds_dir))   mmo <- add_sim_matrix(mmo, m2ds_dir,   "m2ds.sim")

  if (is.null(cos_dir) && is.null(dreams_dir) && is.null(m2ds_dir)) {
    stop("Please provide at least one valid directory.")
  }

  return(mmo)
}


# DEPRECATED: superseded by AddCustomSim().
# AddCustomDist() expects dissimilarity values (0 = identical, 1 = maximally
# different) and stores them in .dissim slots. AddCustomSim() expects
# similarity values (1 = identical, 0 = maximally different) and stores
# them in .sim slots. Retained until downstream functions are updated.
#' Add a custom feature distance matrix to the \code{mmo}
#'
#' Stores any user-supplied pairwise feature distance matrix in the \code{mmo}
#' so it can be used by \code{GetBetaDiversity()}, \code{GetAlphaDiversity()},
#' and \code{HeatmapPlot()} via the \code{distance} argument.
#'
#' @details
#' \strong{Required matrix format:}
#' \itemize{
#'   \item A square numeric matrix with equal row and column names.
#'   \item Row and column names must be feature IDs matching the \code{id} column
#'         in \code{mmo$feature_data} (e.g. \code{"1234.56_2.34"} or whatever
#'         identifier MZmine assigned).
#'   \item Values must be \strong{dissimilarities} in the range \code{[0, 1]},
#'         where 0 means identical and 1 means maximally different. If your
#'         source data are similarities (e.g. Tanimoto coefficients, cosine
#'         scores), convert first with \code{dist_matrix <- 1 - sim_matrix}.
#'   \item The diagonal should be 0 (a feature is identical to itself).
#'   \item The matrix should be symmetric (\code{dist[i,j] == dist[j,i]}).
#'         Non-symmetric matrices are accepted but may produce unexpected results
#'         in downstream ordination and Hill number calculations.
#'   \item Features not present in \code{mmo$feature_data} are retained in the
#'         matrix but will be silently ignored during analysis. Features in
#'         \code{mmo$feature_data} that are absent from the matrix will be
#'         excluded from structure-aware calculations.
#' }
#'
#' @param mmo The \code{mmo}
#' @param dist_matrix A square numeric matrix of pairwise feature dissimilarities
#'   (see Details for format requirements).
#' @param name A short string used to identify this distance in downstream
#'   functions (e.g. \code{"tanimoto"}, \code{"npc"}). The matrix is stored as
#'   \code{mmo$<name>.dissim} and referenced by passing \code{distance = "<name>"}
#'   to \code{GetBetaDiversity()}, \code{GetAlphaDiversity()}, etc.
#'   Must not conflict with existing names: \code{"dreams"}, \code{"cosine"},
#'   \code{"m2ds"}.
#' @return The \code{mmo} with the new distance matrix stored in
#'   \code{mmo$<name>.dissim}.
#' @export
#' @examplesIf FALSE
#' # Example: add a Tanimoto-based structural distance
#' # tanimoto_sim is a feature x feature similarity matrix from an external tool
#' tanimoto_dist <- 1 - tanimoto_sim
#' mmo <- AddCustomDist(mmo, dist_matrix = tanimoto_dist, name = "tanimoto")
#' # Now use it anywhere a distance argument is accepted:
#' beta <- GetBetaDiversity(mmo, method = "CSCS", distance = "tanimoto")
#' alpha <- GetAlphaDiversity(mmo, mode = "weighted", distance = "tanimoto")
AddCustomDist <- function(mmo, dist_matrix, name) {
  if (!is.matrix(dist_matrix) || !is.numeric(dist_matrix)) {
    stop("dist_matrix must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(dist_matrix) != ncol(dist_matrix)) {
    stop("dist_matrix must be square (nrow == ncol).", call. = FALSE)
  }
  if (is.null(rownames(dist_matrix)) || is.null(colnames(dist_matrix))) {
    stop("dist_matrix must have row and column names matching feature IDs in mmo$feature_data$id.", call. = FALSE)
  }
  if (!identical(rownames(dist_matrix), colnames(dist_matrix))) {
    stop("dist_matrix row names and column names must be identical (same features, same order).", call. = FALSE)
  }
  reserved <- c("dreams", "cosine", "m2ds")
  if (name %in% reserved) {
    warning("'", name, "' is a reserved distance name used by AddChemDist(). ",
            "This will overwrite the existing matrix.", call. = FALSE)
  }
  slot_name <- paste0(name, ".dissim")
  mmo[[slot_name]] <- dist_matrix
  message(slot_name, " added to mmo (", nrow(dist_matrix), " x ", ncol(dist_matrix), " features)")
  return(mmo)
}


#' Add a custom feature similarity matrix to the \code{mmo}
#'
#' Stores any user-supplied pairwise feature similarity matrix in the \code{mmo}
#' so it can be used by \code{GetBetaDiversity()}, \code{GetAlphaDiversity()},
#' and \code{HeatmapPlot()} via the \code{distance} argument. Accepts either
#' a dense numeric matrix or a sparse \code{dgCMatrix} (from the Matrix package).
#'
#' @details
#' \strong{Required matrix format:}
#' \itemize{
#'   \item A square numeric matrix (or sparse \code{dgCMatrix}) with equal row
#'         and column names.
#'   \item Row and column names must be feature IDs matching the \code{id} column
#'         in \code{mmo$feature_data}.
#'   \item Values must be \strong{similarities} in the range \code{[0, 1]},
#'         where 1 means identical and 0 means maximally different.
#'   \item The diagonal should be 1 (a feature is identical to itself).
#'   \item The matrix should be symmetric (\code{sim[i,j] == sim[j,i]}).
#'   \item Features not present in \code{mmo$feature_data} are retained but
#'         silently ignored during analysis.
#' }
#'
#' @param mmo The \code{mmo}
#' @param sim_matrix A square numeric matrix or sparse Matrix of pairwise feature
#'   similarities (see Details for format requirements).
#' @param name A short string used to identify this matrix in downstream
#'   functions (e.g. \code{"tanimoto"}, \code{"npc"}). Stored as
#'   \code{mmo$<name>.sim} and referenced via \code{distance = "<name>"}.
#'   Must not conflict with built-in names: \code{"dreams"}, \code{"cosine"},
#'   \code{"m2ds"}.
#' @return The \code{mmo} with the similarity matrix stored in
#'   \code{mmo$<name>.sim}.
#' @export
#' @examplesIf FALSE
#' # tanimoto_sim is a feature x feature similarity matrix from an external tool
#' mmo <- AddCustomSim(mmo, sim_matrix = tanimoto_sim, name = "tanimoto")
#' # Now use it anywhere a distance argument is accepted:
#' beta <- GetBetaDiversity(mmo, method = "CSCS", distance = "tanimoto")
#' alpha <- GetAlphaDiversity(mmo, mode = "weighted", distance = "tanimoto")
AddCustomSim <- function(mmo, sim_matrix, name) {
  if (!is.matrix(sim_matrix) && !inherits(sim_matrix, "Matrix")) {
    stop("sim_matrix must be a numeric matrix or sparse Matrix (e.g. dgCMatrix).",
         call. = FALSE)
  }
  if (is.matrix(sim_matrix) && !is.numeric(sim_matrix)) {
    stop("sim_matrix must contain numeric values.", call. = FALSE)
  }
  if (nrow(sim_matrix) != ncol(sim_matrix)) {
    stop("sim_matrix must be square (nrow == ncol).", call. = FALSE)
  }
  if (is.null(rownames(sim_matrix)) || is.null(colnames(sim_matrix))) {
    stop(
      "sim_matrix must have row and column names matching feature IDs ",
      "in mmo$feature_data$id.",
      call. = FALSE
    )
  }
  if (!identical(rownames(sim_matrix), colnames(sim_matrix))) {
    stop(
      "sim_matrix row names and column names must be identical ",
      "(same features, same order).",
      call. = FALSE
    )
  }
  reserved <- c("dreams", "cosine", "m2ds")
  if (name %in% reserved) {
    warning("'", name, "' is a reserved name used by AddChemSim(). ",
            "This will overwrite the existing matrix.", call. = FALSE)
  }
  slot_name <- paste0(name, ".sim")
  mmo[[slot_name]] <- sim_matrix
  message(slot_name, " added to mmo (",
          nrow(sim_matrix), " x ", ncol(sim_matrix), " features)")
  return(mmo)
}


#' Reorder samples in the \code{mmo} based on group order
#'
#' This function reorders the samples in the \code{mmo} based on a specified group order.
#' The function updates the order of samples in the feature data, log-normalized data, z-score data, and mean-centered data.
#' Use this function before plotting heatmaps or other visualizations to ensure consistent group ordering.
#'
#' @param mmo The \code{mmo}
#' @param group_order A vector specifying the desired order of groups
#' @return The \code{mmo} with reordered samples
#' @export
#' @examplesIf FALSE
#' mmo <- ReorderGroups(mmo, group_order = c("Control", "Treatment1", "Treatment2"))
ReorderGroups <- function(mmo, group_order) {
  metadata <- mmo$metadata
  # Get sample names in the specified group order
  ordered_samples <- unlist(lapply(group_order, function(g) metadata$sample[metadata$group == g]))
  # Reorder columns: id, feature, then ordered samples
  mmo$feature_data <- mmo$feature_data |>
    dplyr::select(.data$id, .data$feature, all_of(ordered_samples))
  mmo$log <- mmo$log |>
    dplyr::select(.data$id, .data$feature, all_of(ordered_samples))
  mmo$zscore <- mmo$zscore |>
    dplyr::select(.data$id, .data$feature, all_of(ordered_samples))
  mmo$meancentered <- mmo$meancentered |>
    dplyr::select(.data$id, .data$feature, all_of(ordered_samples))
  return(mmo)
}

#' Filter an MGF file to keep only spectra for features present in \code{mmo$feature_data$id}
#'
#' Create a new \code{.mgf} file that contains only spectra (ION blocks) whose
#' \code{FEATURE_ID} occurs in \code{mmo$feature_data$id}. This is useful for keeping
#' your spectral library in sync with the features currently stored in a \code{mmo} 
#' (e.g., after subsetting, filtering, or rebuilding a \code{mmo}).
#'
#' Each spectrum in an MGF is represented by a \code{BEGIN IONS ... END IONS} block.
#' This function keeps or discards entire blocks based on the integer value in the
#' \code{FEATURE_ID=} header line. All blocks (MS1 and MS2) are kept for a retained
#' feature, including multiple MS2 blocks if present.
#'
#' @param mmo An eCOMET \code{mmo} containing a required \code{feature_data}
#'   table with an \code{id} column (\code{mmo$feature_data$id}).
#'
#' @param mgf_path Character. Path to the input \code{.mgf} file.
#'
#' @param output_path Character or NULL. Path to write the filtered \code{.mgf}.
#'   If \code{NULL} (default), the output name is derived from \code{mgf_path} by
#'   appending \code{_filtered} before the \code{.mgf} extension (or adding
#'   \code{_filtered.mgf} if no extension is present).
#'
#' @param chunk_lines Integer. Number of lines read per iteration. Larger values are
#'   typically faster but use more memory. Default is \code{100000L}.
#'
#' @param verbose Logical. If \code{TRUE} (default), prints a short progress summary
#'   and final counts.
#'
#' @return Invisibly returns a list with:
#' \itemize{
#'   \item \code{output_path}: path to the filtered MGF
#'   \item \code{blocks_total}: total \code{BEGIN IONS} blocks encountered
#'   \item \code{blocks_kept}: number of blocks written to \code{output_path}
#'   \item \code{lines_read}: total lines read from \code{mgf_path}
#' }
#'
#' @examples
#' \dontrun{
#' # Write "<input>_filtered.mgf" containing only FEATURE_IDs present in mmo$feature_data$id
#' filter_mgf_to_mmo(mmo, "spectra.mgf")
#'
#' # Write to a custom file name
#' filter_mgf_to_mmo(mmo, "spectra.mgf", output_path = "spectra_mmo_only.mgf")
#' }
#'
#' @export
filter_mgf_to_mmo <- function(
    mmo,
    mgf_path,
    output_path = NULL,
    chunk_lines = 100000L,
    verbose = TRUE
) {
  # --- checks ---
  if (is.null(mmo) || !is.list(mmo)) {
    stop("`mmo` must be a list-like ecomet object.", call. = FALSE)
  }
  if (!("feature_data" %in% names(mmo)) || !is.data.frame(mmo$feature_data)) {
    stop("`mmo$feature_data` must exist and be a data.frame.", call. = FALSE)
  }
  if (!("id" %in% names(mmo$feature_data))) {
    stop("`mmo$feature_data$id` must exist.", call. = FALSE)
  }
  if (!is.character(mgf_path) || length(mgf_path) != 1L || !nzchar(mgf_path) || !file.exists(mgf_path)) {
    stop("`mgf_path` must be an existing file path.", call. = FALSE)
  }
  if (!is.null(output_path)) {
    if (!is.character(output_path) || length(output_path) != 1L || !nzchar(output_path)) {
      stop("`output_path` must be NULL or a non-empty character string.", call. = FALSE)
    }
  } else {
    if (grepl("\\.mgf$", mgf_path, ignore.case = TRUE)) {
      output_path <- sub("\\.mgf$", "_filtered.mgf", mgf_path, ignore.case = TRUE)
    } else {
      output_path <- paste0(mgf_path, "_filtered.mgf")
    }
  }

  # --- build ID lookup (fast membership checks) ---
  ids <- mmo$feature_data$id
  ids_num <- suppressWarnings(as.integer(as.character(ids)))
  ids_num <- ids_num[!is.na(ids_num)]
  if (length(ids_num) == 0) {
    stop(
      "No usable (integer) IDs found in `mmo$feature_data$id`.\n",
      "MGF filtering expects FEATURE_ID-like integers.",
      call. = FALSE
    )
  }
  keep_ids <- unique(ids_num)

  # Environment hash set for O(1) membership
  keep_env <- new.env(parent = emptyenv(), hash = TRUE)
  for (k in keep_ids) keep_env[[as.character(k)]] <- TRUE

  # --- streaming filter ---
  in_con  <- file(mgf_path, open = "r", encoding = "UTF-8")
  on.exit(try(close(in_con), silent = TRUE), add = TRUE)
  out_con <- file(output_path, open = "w", encoding = "UTF-8")
  on.exit(try(close(out_con), silent = TRUE), add = TRUE)

  begin_line <- "BEGIN IONS"
  end_line   <- "END IONS"
  key_prefix <- "FEATURE_ID="

  inside_block <- FALSE
  block_lines <- character(0)
  decided_keep <- NA  # logical once FEATURE_ID parsed
  blocks_total <- 0L
  blocks_kept  <- 0L
  lines_read_total <- 0L

  parse_feature_id <- function(line) {
    if (!startsWith(line, key_prefix)) return(NA_integer_)
    val <- substr(line, nchar(key_prefix) + 1L, nchar(line))
    suppressWarnings(as.integer(val))
  }

  if (verbose) {
    message(
      "filter_mgf_to_mmo():\n",
      "- Input MGF:  ", mgf_path, "\n",
      "- Output MGF: ", output_path, "\n",
      "- Using IDs:  mmo$feature_data$id (n unique = ", length(keep_ids), ")\n",
      "- Chunk lines: ", chunk_lines
    )
  }

  repeat {
    chunk <- readLines(in_con, n = chunk_lines, warn = FALSE)
    if (length(chunk) == 0) break
    lines_read_total <- lines_read_total + length(chunk)

    for (line in chunk) {
      if (!inside_block) {
        if (identical(line, begin_line)) {
          inside_block <- TRUE
          block_lines <- line
          decided_keep <- NA
        }
        next
      }

      block_lines <- c(block_lines, line)

      if (is.na(decided_keep)) {
        fid <- parse_feature_id(line)
        if (!is.na(fid)) {
          decided_keep <- isTRUE(keep_env[[as.character(fid)]])
        }
      }

      if (identical(line, end_line)) {
        blocks_total <- blocks_total + 1L

        if (isTRUE(decided_keep)) {
          writeLines(block_lines, out_con, sep = "\n")
          blocks_kept <- blocks_kept + 1L
        }

        inside_block <- FALSE
        block_lines <- character(0)
        decided_keep <- NA
      }
    }

    if (verbose && blocks_total > 0L && (blocks_total %% 5000L == 0L)) {
      message("... processed blocks: ", blocks_total, " | kept: ", blocks_kept)
    }
  }

  if (verbose) {
    message(
      "Done.\n",
      "- Lines read:   ", format(lines_read_total, big.mark = ","), "\n",
      "- Blocks total: ", format(blocks_total, big.mark = ","), "\n",
      "- Blocks kept:  ", format(blocks_kept, big.mark = ","), "\n",
      "- Kept fraction: ", if (blocks_total > 0) sprintf("%.3f", blocks_kept / blocks_total) else "NA"
    )
  }

  invisible(list(
    output_path = output_path,
    blocks_total = blocks_total,
    blocks_kept = blocks_kept,
    lines_read = lines_read_total
  ))
}



#' Annotate \code{mmo$feature_info} with MS2 presence and MS2 block counts from an MGF
#'
#' Scan an \code{.mgf} file and summarize MS/MS availability for each feature in
#' \code{mmo$feature_info}. The function adds two columns:
#' \itemize{
#'   \item \code{ms2}: \code{TRUE} if the MGF contains at least one \code{MSLEVEL=2}
#'         block for that \code{id}; otherwise \code{FALSE}.
#'   \item \code{count_ms2}: number of \code{MSLEVEL=2} blocks for that \code{id}.
#' }
#'
#' This is useful for quickly identifying which features have MS/MS spectra available
#' (and how many replicate MS2 spectra exist) before downstream annotation, networking,
#' or library-building steps.
#'
#' @param mmo An eCOMET \code{mmo} containing a required \code{feature_info}
#'   table with an \code{id} column (\code{mmo$feature_info$id}).
#'
#' @param mgf_path Character. Path to the input \code{.mgf} file.
#'
#' @param chunk_lines Integer. Number of lines read per iteration. Larger values are
#'   typically faster but use more memory. Default is \code{100000L}.
#'
#' @param overwrite Logical. If \code{FALSE} (default) and \code{ms2} and/or \code{count_ms2}
#'   already exist in \code{mmo$feature_info}, the function errors. Set \code{overwrite = TRUE}
#'   to replace existing columns.
#'
#' @param verbose Logical. If \code{TRUE} (default), prints a brief summary of how many
#'   MS2 blocks were found and how many features have MS2.
#'
#' @return The updated \code{mmo} with \code{mmo$feature_info$ms2} and
#'   \code{mmo$feature_info$count_ms2} added (or overwritten if \code{overwrite = TRUE}).
#'
#' @examples
#' \dontrun{
#' # Add ms2 + count_ms2 columns to mmo$feature_info
#' mmo <- annotate_feature_info_ms2_from_mgf(mmo,
#'        "spectra.mgf")
#'
#' # Overwrite existing columns if you re-run on a different MGF
#' mmo <- annotate_feature_info_ms2_from_mgf(mmo,
#'        "spectra_new.mgf", overwrite = TRUE)
#' }
#'
#' @export
annotate_feature_info_ms2_from_mgf <- function(
    mmo,
    mgf_path,
    chunk_lines = 100000L,
    overwrite = FALSE,
    verbose = TRUE
) {
  # --- checks ---
  if (is.null(mmo) || !is.list(mmo)) {
    stop("`mmo` must be a list-like ecomet object.", call. = FALSE)
  }
  if (!("feature_info" %in% names(mmo)) || !is.data.frame(mmo$feature_info)) {
    stop("`mmo$feature_info` must exist and be a data.frame.", call. = FALSE)
  }
  if (!("id" %in% names(mmo$feature_info))) {
    stop("`mmo$feature_info$id` must exist.", call. = FALSE)
  }
  if (!is.character(mgf_path) || length(mgf_path) != 1L || !nzchar(mgf_path) || !file.exists(mgf_path)) {
    stop("`mgf_path` must be an existing file path.", call. = FALSE)
  }
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("`overwrite` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!overwrite) {
    if ("ms2" %in% names(mmo$feature_info) || "count_ms2" %in% names(mmo$feature_info)) {
      stop(
        "Columns already exist in mmo$feature_info: ",
        paste(intersect(c("ms2", "count_ms2"), names(mmo$feature_info)), collapse = ", "),
        "\nSet overwrite = TRUE to overwrite them.",
        call. = FALSE
      )
    }
  }

  # IDs we care about
  ids <- mmo$feature_info$id
  ids_num <- suppressWarnings(as.integer(as.character(ids)))
  ok <- !is.na(ids_num)
  if (!any(ok)) {
    stop("No usable (integer) IDs found in `mmo$feature_info$id`.", call. = FALSE)
  }
  ids_num <- ids_num[ok]
  keep_ids <- unique(ids_num)

  # Environment for fast membership: only count IDs present in feature_info
  id_env <- new.env(parent = emptyenv(), hash = TRUE)
  for (k in keep_ids) id_env[[as.character(k)]] <- TRUE

  # Counts stored in an environment too (sparse)
  count_env <- new.env(parent = emptyenv(), hash = TRUE)

  # --- streaming parse ---
  in_con <- file(mgf_path, open = "r", encoding = "UTF-8")
  on.exit(try(close(in_con), silent = TRUE), add = TRUE)

  begin_line <- "BEGIN IONS"
  end_line   <- "END IONS"
  fid_prefix <- "FEATURE_ID="
  ms_prefix  <- "MSLEVEL="

  inside_block <- FALSE
  block_fid <- NA_integer_
  block_is_ms2 <- FALSE
  saw_fid <- FALSE
  saw_ms <- FALSE

  blocks_total <- 0L
  blocks_ms2_total <- 0L
  lines_read_total <- 0L

  parse_int_after_prefix <- function(line, prefix) {
    if (!startsWith(line, prefix)) return(NA_integer_)
    val <- substr(line, nchar(prefix) + 1L, nchar(line))
    suppressWarnings(as.integer(val))
  }

  if (verbose) {
    message(
      "annotate_feature_info_ms2_from_mgf():\n",
      "- MGF: ", mgf_path, "\n",
      "- Feature IDs: mmo$feature_info$id (n unique = ", length(keep_ids), ")\n",
      "- Chunk lines: ", chunk_lines
    )
  }

  repeat {
    chunk <- readLines(in_con, n = chunk_lines, warn = FALSE)
    if (length(chunk) == 0) break
    lines_read_total <- lines_read_total + length(chunk)

    for (line in chunk) {
      if (!inside_block) {
        if (identical(line, begin_line)) {
          inside_block <- TRUE
          block_fid <- NA_integer_
          block_is_ms2 <- FALSE
          saw_fid <- FALSE
          saw_ms <- FALSE
        }
        next
      }

      # Inside a block: look only for FEATURE_ID and MSLEVEL
      if (!saw_fid && startsWith(line, fid_prefix)) {
        fid <- parse_int_after_prefix(line, fid_prefix)
        if (!is.na(fid)) {
          block_fid <- fid
          saw_fid <- TRUE
        }
        next
      }

      if (!saw_ms && startsWith(line, ms_prefix)) {
        ms <- parse_int_after_prefix(line, ms_prefix)
        if (!is.na(ms) && ms == 2L) {
          block_is_ms2 <- TRUE
        }
        saw_ms <- TRUE
        next
      }

      if (identical(line, end_line)) {
        blocks_total <- blocks_total + 1L

        # Count only if:
        # - it's MS2
        # - FEATURE_ID parsed
        # - FEATURE_ID is in our feature_info set
        if (block_is_ms2 && !is.na(block_fid) && isTRUE(id_env[[as.character(block_fid)]])) {
          blocks_ms2_total <- blocks_ms2_total + 1L
          key <- as.character(block_fid)
          cur <- count_env[[key]]
          if (is.null(cur)) cur <- 0L
          count_env[[key]] <- cur + 1L
        }

        inside_block <- FALSE
      }
    }
  }

  # Build count vector aligned to feature_info rows
  fi <- mmo$feature_info
  fi_ids_num <- suppressWarnings(as.integer(as.character(fi$id)))

  count_ms2 <- integer(nrow(fi))
  for (i in seq_len(nrow(fi))) {
    fid <- fi_ids_num[[i]]
    if (!is.na(fid)) {
      v <- count_env[[as.character(fid)]]
      if (!is.null(v)) count_ms2[[i]] <- as.integer(v)
    }
  }
  ms2 <- count_ms2 > 0L

  fi$ms2 <- ms2
  fi$count_ms2 <- count_ms2
  mmo$feature_info <- fi

  if (verbose) {
    message(
      "Done.\n",
      "- Lines read: ", format(lines_read_total, big.mark = ","), "\n",
      "- Blocks total: ", format(blocks_total, big.mark = ","), "\n",
      "- MS2 blocks counted (in feature_info): ", format(blocks_ms2_total, big.mark = ","), "\n",
      "- Features with MS2: ", sum(ms2, na.rm = TRUE), " / ", nrow(fi)
    )
  }

  mmo
}






#' Filter a \code{mmo} by samples, groups, and/or features
#'
#' @description
#' Subset all components of a \code{mmo} (feature tables, metadata, distance
#' matrices, and annotations) to a given set of samples, groups, and/or
#' feature IDs. Filtering is applied consistently across all slots present
#' in \code{mmo}.
#'
#' Optionally, if an \code{mgf_path} is provided, a filtered MGF is written
#' containing only spectra for the retained features (using \code{filter_mgf_to_mmo()}).
#'
#' @param mmo A list-like \code{mmo} as returned by \code{GetMZmineFeature()}.
#' @param sample_list Optional character vector of sample IDs (matching
#'   \code{sample_col} in \code{mmo$metadata}) to retain.
#' @param group_list Optional character vector of group labels (matching
#'   \code{group_col} in \code{mmo$metadata}) to retain. Mutually exclusive with
#'   \code{sample_list}.
#' @param id_list Optional character vector of feature IDs to retain.
#'   If \code{NULL}, features are determined from \code{feature_data} and optionally
#'   filtered by \code{drop_empty_feat}.
#' @param sample_col Column name in \code{mmo$metadata} containing sample IDs.
#'   Default is \code{"sample"}.
#' @param group_col Column name in \code{mmo$metadata} containing group labels.
#'   Default is \code{"group"}.
#' @param drop_empty_feat Logical; if \code{TRUE} (default) drop features with no
#'   non-zero values in the retained samples.
#' @param empty_threshold Optional numeric threshold used to define "empty"
#'   features. If \code{NULL} (default), the smallest positive, non-NA intensity
#'   in the retained samples is used. Features are kept if they have at least
#'   one value > threshold across retained samples.
#'
#' @param mgf_path Optional character. If provided, an MGF file will be filtered
#'   to retained features using \code{filter_mgf_to_mmo()}.
#' @param output_path Character or NULL. Passed to \code{filter_mgf_to_mmo()}.
#'   If \code{NULL} (default), output is \code{"<input>_filtered.mgf"}.
#'
#' @return A filtered \code{mmo} with the same structure as \code{mmo}, but
#'   restricted to the requested samples / groups / features. If \code{mgf_path}
#'   is provided, the returned object also includes \code{mmo_filtered$mgf_filtered_path}.
#' @export
#'
#' @examples
#' \dontrun{
#' mmo_sub <- filter_mmo(mmo, group_list = c("Species1", "Species2"))
#'
#' # Also write a filtered mgf:
#' mmo_sub <- filter_mmo(mmo, group_list = c("Species1", "Species2"),
#'                       mgf_path = "spectra.mgf")
#' }
filter_mmo <- function(mmo,
                       sample_list  = NULL,
                       group_list   = NULL,
                       id_list = NULL,
                       sample_col   = "sample",
                       group_col    = "group",
                       drop_empty_feat   = TRUE,
                       empty_threshold   = NULL,
                       mgf_path          = NULL,
                       output_path       = NULL) {

  # ------------------------------------------------------------------
  # 0. Basic checks on arguments
  # ------------------------------------------------------------------
  if (is.null(sample_list) && is.null(group_list) && is.null(id_list)) {
    stop("Must provide at least one of 'sample_list', 'group_list', or 'id_list'.", call. = FALSE)
  }
  if (!is.null(sample_list) && !is.null(group_list)) {
    stop("Provide exactly one of 'sample_list' or 'group_list', not both.", call. = FALSE)
  }
  if (!"metadata" %in% names(mmo)) stop("mmo$metadata must be present.", call. = FALSE)
  if (!"feature_data" %in% names(mmo)) stop("mmo$feature_data must be present.", call. = FALSE)

  meta <- mmo$metadata
  if (!is.data.frame(meta)) stop("mmo$metadata must be a data.frame.", call. = FALSE)
  if (!sample_col %in% names(meta)) stop(sprintf("sample_col '%s' not present in metadata.", sample_col), call. = FALSE)
  if (!is.null(group_list) && !group_col %in% names(meta)) {
    stop(sprintf("group_col '%s' not present in metadata.", group_col), call. = FALSE)
  }

  # ------------------------------------------------------------------
  # 1. Determine which samples to keep
  # ------------------------------------------------------------------
  all_samples <- as.character(meta[[sample_col]])

  if (!is.null(sample_list)) {
    keep_samples <- intersect(all_samples, as.character(sample_list))
    if (length(keep_samples) == 0L) {
      stop("None of the supplied sample IDs in 'sample_list' are present in the mmo.", call. = FALSE)
    }
  } else if (!is.null(group_list)) {
    idx <- meta[[group_col]] %in% group_list
    keep_samples <- as.character(meta[[sample_col]][idx])
    if (length(keep_samples) == 0L) {
      stop("No samples found for the requested 'group_list'.", call. = FALSE)
    }
  } else {
    keep_samples <- all_samples
  }

  # ------------------------------------------------------------------
  # 2. Filter feature_data and decide which features to keep
  # ------------------------------------------------------------------
  feat_abund <- mmo[["feature_data"]]
  if (!is.data.frame(feat_abund) || !all(c("id", "feature") %in% names(feat_abund))) {
    stop("mmo$feature_data must be a data.frame with columns 'id' and 'feature'.", call. = FALSE)
  }

  sample_cols_present <- intersect(keep_samples, names(feat_abund))
  col_keep <- c("id", "feature", sample_cols_present)
  feat_abund_filtered <- feat_abund[, col_keep, drop = FALSE]

  auto_features <- as.character(feat_abund_filtered$id)

  if (drop_empty_feat && length(sample_cols_present) > 0L) {
    num_mat <- as.matrix(
      suppressWarnings(
        sapply(feat_abund_filtered[, sample_cols_present, drop = FALSE], as.numeric)
      )
    )
    if (!is.matrix(num_mat)) num_mat <- matrix(num_mat, nrow = nrow(feat_abund_filtered))

    positive_vals <- num_mat[num_mat > 0 & !is.na(num_mat)]

    if (is.null(empty_threshold)) {
      if (length(positive_vals) == 0L) {
        threshold <- Inf
        nonempty_idx <- rep(FALSE, nrow(num_mat))
        message("drop_empty_feat = TRUE, but no positive intensities found; all features are treated as empty.")
      } else {
        min_val <- min(positive_vals)
        threshold <- min_val
        message(sprintf("Removing empty features based on min_val = %g", min_val))
        nonempty_idx <- rowSums(num_mat > threshold, na.rm = TRUE) > 0
      }
    } else {
      threshold <- empty_threshold
      message(sprintf("Removing empty features using user threshold = %g", threshold))
      nonempty_idx <- rowSums(num_mat > threshold, na.rm = TRUE) > 0
    }

    auto_features <- auto_features[nonempty_idx]
    feat_abund_filtered <- feat_abund_filtered[nonempty_idx, , drop = FALSE]
  }

  if (!is.null(id_list)) {
    id_list <- as.character(id_list)
    features_keep <- intersect(auto_features, id_list)
    if (length(features_keep) == 0L) {
      stop("No overlap between 'id_list' and features present in mmo$feature_data.", call. = FALSE)
    }
    feat_abund_filtered <- feat_abund_filtered[feat_abund_filtered$id %in% features_keep, , drop = FALSE]
  } else {
    features_keep <- auto_features
  }

  # ------------------------------------------------------------------
  # 3. Helpers to filter common slot types
  # ------------------------------------------------------------------
  filter_feature_sample_table <- function(df, keep_ids, keep_samples) {
    if (!is.data.frame(df)) return(df)
    if (!all(c("id", "feature") %in% names(df))) return(df)
    samp <- intersect(keep_samples, names(df))
    cols <- c("id", "feature", samp)
    out <- df[, cols, drop = FALSE]
    out[out$id %in% keep_ids, , drop = FALSE]
  }

  filter_feature_df_by_id <- function(df, keep_ids) {
    if (!is.data.frame(df)) return(df)
    if (!("id" %in% names(df))) return(df)
    df[df$id %in% keep_ids, , drop = FALSE]
  }

  filter_square_feature_matrix <- function(mat, keep_ids) {
    if (!is.matrix(mat) && !inherits(mat, "Matrix")) return(mat)
    rn <- rownames(mat); cn <- colnames(mat)
    if (is.null(rn) || is.null(cn)) return(mat)
    # only treat as square feature-feature matrix if row/col names overlap a lot with keep_ids
    if (nrow(mat) != ncol(mat)) return(mat)
    keep_r <- rn %in% keep_ids
    keep_c <- cn %in% keep_ids
    if (!any(keep_r) || !any(keep_c)) return(mat)
    mat[keep_r, keep_c, drop = FALSE]
  }

  # ------------------------------------------------------------------
  # 4. Build filtered mmo (preserve *all* slots, filtered where appropriate)
  # ------------------------------------------------------------------
  mmo_filtered <- mmo

  # Required slots
  mmo_filtered$feature_data <- feat_abund_filtered
  mmo_filtered$metadata <- meta[meta[[sample_col]] %in% keep_samples, , drop = FALSE]

  # Standard featurexsample tables if present
  for (nm in c("log", "zscore", "meancentered")) {
    if (nm %in% names(mmo_filtered)) {
      mmo_filtered[[nm]] <- filter_feature_sample_table(mmo_filtered[[nm]], features_keep, keep_samples)
    }
  }

  # Filter feature_info if present
  if ("feature_info" %in% names(mmo_filtered)) {
    mmo_filtered$feature_info <- filter_feature_df_by_id(mmo_filtered$feature_info, features_keep)
  }

  # Pairwise tables (feature-level data.frames with id)
  if ("pairwise" %in% names(mmo_filtered)) {
    mmo_filtered$pairwise <- filter_feature_df_by_id(mmo_filtered$pairwise, features_keep)
  }

  # Any "annotation-like" feature-level data.frame with an id column:
  # - this catches sirius_annot AND any sirius_annot_filtered_* you've created,
  #   plus other feature-level tables you add in the future.
  for (nm in names(mmo_filtered)) {
    obj <- mmo_filtered[[nm]]

    # Skip those we handled explicitly above
    if (nm %in% c("feature_data", "metadata", "log", "zscore", "meancentered", "feature_info", "pairwise")) next

    # Feature-level data.frames: subset rows by id
    if (is.data.frame(obj) && "id" %in% names(obj)) {
      mmo_filtered[[nm]] <- obj[obj$id %in% features_keep, , drop = FALSE]
      next
    }

    # Square feature-feature matrices: subset by row/col names
    if (is.matrix(obj) || inherits(obj, "Matrix")) {
      mmo_filtered[[nm]] <- filter_square_feature_matrix(obj, features_keep)
      next
    }
  }

  # Preserve class
  class(mmo_filtered) <- class(mmo)

  # ------------------------------------------------------------------
  # 5. Optional: filter an MGF in the same call
  # ------------------------------------------------------------------
  if (!is.null(mgf_path)) {
    if (!is.character(mgf_path) || length(mgf_path) != 1L || !nzchar(mgf_path) || !file.exists(mgf_path)) {
      stop("`mgf_path` must be a single existing file path.", call. = FALSE)
    }
    mgf_res <- filter_mgf_to_mmo(
      mmo = mmo_filtered,
      mgf_path = mgf_path,
      output_path = output_path,
      verbose = FALSE
    )
    mmo_filtered$mgf_filtered_path <- mgf_res$output_path
  }

  message("mmo was subset")
  message(paste0("Feature number: ", nrow(mmo_filtered$feature_data)))
  message(paste0(nrow(mmo_filtered$metadata), " samples in ",
                 length(unique(mmo_filtered$metadata[[group_col]])), " groups"))

  mmo_filtered
}

### Pool Rows by Group

#' pool_mmo_by_group
#'
#' Pools sample columns within each group into one pseudo-sample per group.
#' Feature rows are preserved (no filtering of features).
#' Keeps all other slots of the \code{mmo} unchanged by copying the \code{mmo} first.
#'
#' @param mmo A \code{mmo}
#' @param group_col column in \code{mmo$metadata} used for grouping (default: "group")
#' @return \code{mmo} with feature_data containing one column per group and metadata updated accordingly
#' @export
pool_mmo_by_group <- function(mmo, group_col = "group") {
  if (is.null(mmo$feature_data) || is.null(mmo$metadata)) {
    stop("mmo must contain $feature_data and $metadata.")
  }

  fd <- mmo$feature_data
  md <- mmo$metadata

  if (!("sample" %in% names(md))) stop("mmo$metadata must contain 'sample'.")
  if (!(group_col %in% names(md))) stop("group_col not found in metadata: ", group_col)

  # aligned sample columns (assumes first 2 columns are id + feature)
  sample_cols <- colnames(fd)[-c(1, 2)]
  sample_cols <- intersect(sample_cols, md$sample)
  if (length(sample_cols) < 1) stop("No aligned sample columns to pool.")

  group_map <- as.character(md[[group_col]][match(sample_cols, md$sample)])
  group_map[is.na(group_map) | group_map == ""] <- "ungrouped"

  # stable group order (appearance order in metadata, not alphabetical)
  groups <- unique(group_map)

  x <- as.matrix(fd[, sample_cols, drop = FALSE])

  pooled_mat <- matrix(NA_real_, nrow = nrow(x), ncol = length(groups))
  colnames(pooled_mat) <- groups

  for (i in seq_along(groups)) {
    g <- groups[i]
    cols_g <- sample_cols[group_map == g]
    pooled_mat[, i] <- rowSums(x[, cols_g, drop = FALSE], na.rm = TRUE)
  }

  pooled_fd <- data.frame(
    fd[, 1:2, drop = FALSE],
    as.data.frame(pooled_mat, check.names = FALSE),
    check.names = FALSE
  )
  colnames(pooled_fd)[1:2] <- colnames(fd)[1:2]

  pooled_md <- data.frame(
    sample = groups,
    group = groups,
    stringsAsFactors = FALSE
  )
  if (group_col != "group") pooled_md[[group_col]] <- groups

  # copy whole object, then overwrite only what must change
  out <- mmo
  out$feature_data <- pooled_fd
  out$metadata <- pooled_md
  out
}





########################################################################################
# Define functions for supporting analysis
########################################################################################

#' Retrieve feature data from the \code{mmo}, with normalization options
#'
#' This function retrieves the feature data from the \code{mmo} based on the specified normalization method.
#'
#' @param mmo The \code{mmo}
#' @param normalization The normalization method to use. Options are 'None','PA', 'Log', 'Meancentered', 'Z', or 'Imputed'
#' @return The feature data corresponding to the specified normalization method
#' @export
#' @examplesIf FALSE
#' feature_data <- GetNormFeature(mmo, normalization = 'Log')
GetNormFeature <- function(mmo, normalization){
  if (normalization == 'None'){
    feature <- mmo$feature_data
  } else if (normalization == 'PA'){
    feature <- mmo$feature_presence
  } else if (normalization == 'Log'){
    feature <- mmo$log
  } else if (normalization == 'Meancentered'){
    feature <- mmo$meancentered
  } else if (normalization == 'Z'){
    feature <- mmo$zscore
  } else if (normalization == 'Imputed'){
    feature <- mmo$imputed_feature_data
  } else {
    print('The normalization should be None, PA, Log, Meancentered, Z, or Imputed')
  }
  return(feature)
}


# DEPRECATED: superseded by GetSimMat().
# GetDistanceMat() looks up .dissim slots (created by AddChemDist() /
# AddCustomDist()). GetSimMat() looks up .sim slots (created by AddChemSim()
# / AddCustomSim()). Retained until all downstream callers are updated.
#' Get the distance matrix from the \code{mmo} based on the specified distance metric
#'
#' Retrieve a feature distance matrix from the \code{mmo}
#'
#' Looks up a dissimilarity matrix stored in the \code{mmo} by name. Works with
#' the three built-in matrices added by \code{AddChemDist()} as well as any
#' custom matrix added via \code{AddCustomDist()}.
#'
#' @param mmo The \code{mmo}
#' @param distance Name of the distance matrix to retrieve. Built-in options are
#'   \code{'dreams'}, \code{'cosine'}, and \code{'m2ds'}. Any name passed to
#'   \code{AddCustomDist(mmo, name = ...)} is also valid.
#' @return The distance matrix (numeric matrix with feature ID row/col names).
#' @export
#' @examplesIf FALSE
#' distance_matrix <- GetDistanceMat(mmo, distance = 'dreams')
#' distance_matrix <- GetDistanceMat(mmo, distance = 'tanimoto')  # custom
GetDistanceMat <- function(mmo, distance = 'dreams'){
  slot_name <- paste0(distance, ".dissim")
  distance_matrix <- mmo[[slot_name]]
  if (is.null(distance_matrix)) {
    stop(
      "No distance matrix found for '", distance, "' (looked for mmo$", slot_name, "). ",
      "Use AddChemDist() for built-in types or AddCustomDist() for custom matrices.",
      call. = FALSE
    )
  }
  return(distance_matrix)
}


#' Retrieve a feature similarity matrix from the \code{mmo}
#'
#' Looks up a similarity matrix stored in the \code{mmo} by name. Works with
#' the three built-in matrices added by \code{AddChemSim()} as well as any
#' custom matrix added via \code{AddCustomSim()}. Returns the matrix as-is
#' (dense numeric or sparse dgCMatrix); downstream functions are responsible
#' for converting to distance when needed.
#'
#' @param mmo The \code{mmo}
#' @param distance Name of the similarity matrix to retrieve. Built-in options
#'   are \code{'dreams'}, \code{'cosine'}, and \code{'m2ds'}. Any name passed
#'   to \code{AddCustomSim(mmo, name = ...)} is also valid.
#' @return The similarity matrix (numeric matrix or sparse dgCMatrix with
#'   feature ID row/col names).
#' @export
#' @examplesIf FALSE
#' sim_matrix <- GetSimMat(mmo, distance = 'dreams')
#' sim_matrix <- GetSimMat(mmo, distance = 'tanimoto')  # custom
GetSimMat <- function(mmo, distance = 'dreams') {
  slot_name <- paste0(distance, ".sim")
  sim_matrix <- mmo[[slot_name]]
  if (is.null(sim_matrix)) {
    stop(
      "No similarity matrix found for '", distance,
      "' (looked for mmo$", slot_name, "). ",
      "Use AddChemSim() for built-in types or AddCustomSim() for custom matrices.",
      call. = FALSE
    )
  }
  return(sim_matrix)
}

#' Convert feature names to IDs in the \code{mmo}
#'
#' This function converts feature names to their corresponding IDs in the \code{mmo}.
#'
#' @param mmo The \code{mmo}
#' @param feature_names A vector of feature names to convert
#' @return A vector of feature IDs corresponding to the input feature names
#' @export
#' @examplesIf FALSE
#' feature_ids <- FeatureToID(mmo, feature_names = c("100.0_5.0", "150.0_10.0"))
#' feature_ids <- FeatureToID(mmo, feature_names = mmo$feature_data$feature[1:10])
#' feature_ids <- FeatureToID(mmo,
#'  feature_names = Glucosinolates
#' ) # if Glucosinolates is a vector of feature names
FeatureToID <- function(mmo, feature_names) {
  feature_data <- mmo$feature_data
  feature_ids <- feature_data |>
    filter(.data$feature %in% feature_names) |>
    dplyr::select(.data$feature, .data$id)
  # Match the order of feature_names
  feature_ids <- feature_ids$id[match(feature_names, feature_ids$feature)]
  return(feature_ids)
}

#' Convert feature IDs to names in the \code{mmo}
#'
#' This function converts feature IDs to their corresponding names in the \code{mmo}.
#'
#' @param mmo The \code{mmo}
#' @param feature_ids A vector of feature IDs to convert
#' @return A vector of feature names corresponding to the input feature IDs
#' @export
#' @examplesIf FALSE
#' feature_names <- IDToFeature(mmo, feature_ids = c("1219", "2250", "3360"))
#' feature_names <- IDToFeature(mmo, feature_ids = mmo$feature_data$id[1:10])
#' feature_names <- IDToFeature(mmo,
#'  feature_ids = FeatureToID(mmo, feature_names = Glucosinolates)
#' ) # if Glucosinolates is a vector of feature names
IDToFeature <- function(mmo, feature_ids) {
  feature_data <- mmo$feature_data
  feature_names <- feature_data |>
    filter(.data$id %in% feature_ids) |>
    dplyr::select(.data$id, .data$feature)
  # Match the order of feature_ids
  feature_names <- feature_names$feature[match(feature_ids, feature_names$id)]
  return(feature_names)
}

#' Calculate group means from the \code{mmo}
#'
#' This function calculates and returns a dataframe of mean feature values for each group in the \code{mmo}, with options for normalization and filtering.
#' Use SwitchGroup() to change the grouping variable before running this function.
#'
#' @param mmo The \code{mmo}
#' @param normalization The normalization method to use. Options are 'None', 'Log', 'Meancentered', or 'Z'
#' @param filter_id Boolean to filter features based on a provided list (default: FALSE)
#' @param id_list A vector of feature names to filter (default: NULL)
#' @param filter_group Boolean to filter groups based on a provided list (default: FALSE)
#' @param group_list A vector of group names to filter (default: NULL)
#' @return A data frame containing the mean feature values for each group
#' @export
#' @examplesIf FALSE
#' group_means <- GetGroupMeans(mmo, normalization = 'Log')
#' group_means <- GetGroupMeans(mmo,
#'  normalization = 'None',
#'  filter_id = TRUE, id_list = Glucosinolates
#' ) # if Glucosinolates is a vector of feature names
#' group_means <- GetGroupMeans(mmo,
#'  normalization = 'Z',
#'  filter_group = TRUE,
#'  group_list = c("Control", "Treatment1")
#' )
#' group_means <- GetGroupMeans(mmo,
#'  normalization = 'Meancentered',
#'  filter_id = TRUE, id_list = Glucosinolates,
#'  filter_group = TRUE, group_list = c("Control", "Treatment1")
#' ) # if Glucosinolates is a vector of feature names
GetGroupMeans <- function(mmo, normalization = 'None', filter_id = FALSE, id_list = NULL, filter_group = FALSE, group_list = NULL) {
  if (filter_id || filter_group) {
    mmo <- filter_mmo(mmo, id_list = id_list, group_list = group_list)
  }
  feature_data <- GetNormFeature(mmo, normalization = normalization)
  metadata     <- mmo$metadata

  # Extract numeric matrix (features x samples) — avoids 49M-row pivot_longer
  mat      <- as.matrix(feature_data[, -(1:2), drop = FALSE])
  id_col   <- feature_data$id

  # Only groups that actually have samples, in factor-level order. Matches the
  # columns the previous group_by()/pivot_wider() implementation produced:
  # dplyr drops unused factor levels, so an empty level must not become a column.
  present  <- unique(metadata$group[metadata$sample %in% colnames(mat)])
  groups   <- if (is.factor(metadata$group)) {
    levels(metadata$group)[levels(metadata$group) %in% as.character(present)]
  } else {
    sort(unique(as.character(present)))
  }

  means_list <- lapply(groups, function(g) {
    samps <- metadata$sample[as.character(metadata$group) == g]
    samps <- intersect(samps, colnames(mat))
    if (length(samps) == 0L) return(rep(NA_real_, nrow(mat)))
    rowMeans(mat[, samps, drop = FALSE], na.rm = TRUE)
  })
  names(means_list) <- groups

  result <- data.frame(id = id_col, as.data.frame(means_list, check.names = FALSE),
                       stringsAsFactors = FALSE, check.names = FALSE)
  # Rows are returned sorted by feature id, as group_by(group, id) did, so that
  # heatmap row order and exported tables stay identical to earlier releases.
  result <- result[order(result$id), , drop = FALSE]
  rownames(result) <- NULL
  return(result)
}

#' Calculate log2 fold change for a given control group
#'
#' This function calculates and returns a dataframe of log2 fold change values for each group compared to a specified control group.
#' Takes inputs from GetGroupMeans() function.
#'
#' @param group_means A data frame containing the mean feature values for each group
#' @param control_group The name of the control group to compare against
#' @return A data frame with log2 fold change values for each group compared to the control group
#' @export
#' @examplesIf FALSE
#' fold_change <- GetLog2FoldChange(GetGroupMeans(mmo,
#'                  normalization = 'Log'),
#'                  control_group = 'Control')
GetLog2FoldChange <- function(group_means, control_group) {
  control_means <- group_means[[control_group]]
  fold_change <- group_means |>
    mutate(across(-.data$id, ~ log2(. / control_means)))

  return(fold_change)
}


#' Perform ANOVA and Tukey's HSD test on the \code{mmo}
#'
#' This function performs ANOVA and Tukey's HSD test on the feature data of the \code{mmo},
#' Returns a list of ANOVA results, Tukey's HSD results, Tukey's significance letters, and Dunnett's test results.
#'
#' @param df The data frame containing the feature data and metadata
#' @param formula The formula for the ANOVA test, e.g., "feature ~ group"
#' @return A list containing the ANOVA results, Tukey's HSD results, Tukey's significance letters, and Dunnett's test results
#' @export
#' @examplesIf FALSE
#' anova_results <- anova_tukey_dunnett(df = merged_data, formula = "feature_value ~ group")
anova_tukey_dunnett <- function(df, formula) {
  .require_pkg("DescTools")
  .require_pkg("multcompView")
  aov_res <- aov(as.formula(formula), data = df)
  tukey_res <- TukeyHSD(aov_res)
  tukey_sig <- multcompView::multcompLetters4(aov_res, tukey_res)
  dunnett_res <- DescTools::DunnettTest(as.formula(formula), data = df)
  return(list(aov_res = aov_res, tukey_res = tukey_res, tukey_sig = tukey_sig, dunnett_res = dunnett_res))
}

#' Write results of anova_tukey_dunnett to a CSV file
#'
#' This function writes the results of ANOVA and Tukey's HSD test to a CSV file.
#'
#' @param anova_data A list containing the results of ANOVA and Tukey's HSD test
#' @param outdir The output directory where the results will be saved
#' @param way The type of ANOVA test to perform. Options are 'oneway' or 'twoway'
#' @export
#' @examplesIf FALSE
#' write_anova(anova_data = anova_results, outdir = "anova_tukey_results.csv", way = 'oneway')
#'
write_anova <- function(anova_data, outdir, way='oneway'){
  way_num <- switch(way, oneway = 1, twoway = 3)
  # Perform ANOVA and Tukey HSD
  aov_res <- anova_data$aov_res
  tukey_res <- anova_data$tukey_res
  tukey_sig <- anova_data$tukey_sig[way_num]
  dunnett_res <- anova_data$dunnett_res

  # Save ANOVA and Tukey HSD results
  anova_df <- as.data.frame(summary(aov_res)[[1]])
  anova_df$Comparison <- rownames(anova_df)
  tukey_df <- as.data.frame(tukey_res[way_num])
  tukey_df$Comparison <- rownames(tukey_df)
  sig_letter <- as.data.frame(unlist(tukey_sig))
  sig_letter$Comparison <- rownames(sig_letter)
  dunnett_df <- as.data.frame(dunnett_res[[1]])
  dunnett_df$comp <- rownames(dunnett_df)

  # Create a combined results data frame
  combined_df <- dplyr::bind_rows(
    dplyr::tibble(Test = "ANOVA", anova_df),
    dplyr::tibble(Test = "Tukey", tukey_df),
    dplyr::tibble(Test = 'sig', sig_letter),
    dplyr::tibble(Test = 'Dunnett', dunnett_df)
  )
  readr::write_csv(combined_df, file = outdir)
}

#' Perform PERMANOVA and pairwise comparisons
#'
#' This function performs PERMANOVA on the given data and metadata, with options for filtering groups.
#' It also conducts post-hoc pairwise comparisons and adjusts p-values for multiple testing.
#' The function returns the PERMANOVA results, raw pairwise comparison results, and matrices of adjusted p-values, F values, and R square for pairwise comparisons
#'
#' @param data A data frame or distance matrix for PERMANOVA
#' @param metadata A data frame containing sample metadata, including a 'group' column
#' @param mode The mode of the input data: 'data' for raw data or 'distance' for a distance matrix
#' @param filter_group Boolean to filter groups based on a provided list (default: FALSE)
#' @param group_list A vector of group names to filter (default: NULL)
#' @param permutations The number of permutations for the PERMANOVA test (default: 5000)
#' @return A list containing the PERMANOVA results, raw pairwise comparison results, and matrices of adjusted p-values, F values, and R square for pairwise comparisons
#' @export
#' @examplesIf FALSE
#' permanova_results <- permanova_stat(
#'  data = feature_data, metadata = mmo$metadata,
#'  mode = 'data', filter_group = TRUE, group_list = c("Control", "Treatment1"),
#'  permutations = 5000
#' )
#' permanova_results <- permanova_stat(
#'  data = betadiv, metadata = mmo$metadata,
#'  mode = 'distance', permutations = 10000
#' )
permanova_stat <- function(data, metadata, mode, filter_group = FALSE, group_list = NULL, permutations = 5000){
  .require_pkg("vegan")
  .require_pkg("stringr")
  .require_pkg("pairwiseAdonis")
  #Perform PERMANOVA
  if (mode == 'data'){
    dist_mat <- stats::dist(data)
  } else if (mode == 'distance'){
    dist_mat <- as.dist(data)
  } else {
    stop("Invalid mode for adonis2. Please use 'data' or 'distance'")
  }
  metadata <- metadata[match(rownames(data), metadata$sample), ]
  if (filter_group == TRUE){
    metadata <- metadata |> filter(.data$group %in% group_list)
    dist_mat <- as.dist(as.matrix(dist_mat)[metadata$sample, metadata$sample])
  }
  permanova_res <- vegan::adonis2(dist_mat ~ group, data = metadata, permutations = permutations)
  #Post-hoc
  pairwise_permanova <- pairwiseAdonis::pairwise.adonis2(dist_mat ~ group, data = metadata, permutations = permutations)
  # Multiple test correction for pairwise p-values
  pvals <- sapply(pairwise_permanova[-1], function(x) x[1,5])
  pvals_adj <- p.adjust(pvals, method = "BH")
  # Add adjusted p-values to pairwise_permanova
  for (i in seq_along(pvals_adj)) {
    pairwise_permanova[[i+1]]$padj <- pvals_adj[i]
  }
  #Organize posthoc results for visual
  pairwise_p_matrix <- matrix(NA, nrow = length(unique(metadata$group)), ncol = length(unique(metadata$group)))
  rownames(pairwise_p_matrix) <- unique(metadata$group)
  colnames(pairwise_p_matrix) <- unique(metadata$group)
  for (i in 2:length(pairwise_permanova)){
    group1 <- stringr::str_split(names(pairwise_permanova)[i], '_vs_')[[1]][1]
    group2 <- stringr::str_split(names(pairwise_permanova)[i], '_vs_')[[1]][2]
    pval <- pairwise_permanova[[i]][1,6]
    pairwise_p_matrix[group2, group1] <- pval
  }
  # Organize posthoc F values for visualization
  pairwise_F_matrix <- matrix(NA, nrow = length(unique(metadata$group)), ncol = length(unique(metadata$group)))
  rownames(pairwise_F_matrix) <- unique(metadata$group)
  colnames(pairwise_F_matrix) <- unique(metadata$group)
  for (i in 2:length(pairwise_permanova)){
    group1 <- stringr::str_split(names(pairwise_permanova)[i], '_vs_')[[1]][1]
    group2 <- stringr::str_split(names(pairwise_permanova)[i], '_vs_')[[1]][2]
    Fval <- pairwise_permanova[[i]][1,4]
    pairwise_F_matrix[group2, group1] <- Fval
  }
  # Organize posthoc R^2 values for visualization
  pairwise_R2_matrix <- matrix(NA, nrow = length(unique(metadata$group)), ncol = length(unique(metadata$group)))
  rownames(pairwise_R2_matrix) <- unique(metadata$group)
  colnames(pairwise_R2_matrix) <- unique(metadata$group)
  for (i in 2:length(pairwise_permanova)){
    group1 <- stringr::str_split(names(pairwise_permanova)[i], '_vs_')[[1]][1]
    group2 <- stringr::str_split(names(pairwise_permanova)[i], '_vs_')[[1]][2]
    R2val <- pairwise_permanova[[i]][1,3]
    pairwise_R2_matrix[group2, group1] <- R2val
  }
  return(list(permanova_res = permanova_res, pairwise_raw = pairwise_permanova, pairwise_p_matrix = pairwise_p_matrix, pairwise_F_matrix = pairwise_F_matrix, pairwise_R2_matrix = pairwise_R2_matrix))
}

########################################################################################
# Define functions for pairwise comparison and visualization
########################################################################################

#' Perform pairwise comparison between two groups in the \code{mmo}
#'
#' This function performs pairwise comparison between two groups in the \code{mmo},
#' calculating log2 fold change and adjusted p-values for given comparison of two groups.
#' The function adds the results to the \code{mmo$pairwise} data frame.
#'
#' @param mmo The \code{mmo}
#' @param group1 The name of the nominator group
#' @param group2 The name of the denominator group
#' @param correction The method for multiple comparison correction. Options are 'BH', 'holm', 'bonferroni', etc. Inherits from p.adjust() (default: 'BH')
#' @return The \code{mmo} with pairwise comparison results added to \code{mmo$pairwise}
#' @export
#' @examplesIf FALSE
#' mmo <- PairwiseComp(mmo, group1 = 'Control', group2 = 'Treatment1')
PairwiseComp <- function(mmo, group1, group2, correction = 'BH'){
  if (is.null(mmo$imputed_feature_data)){
    message('No imputed data found. Using feature_data for pairwise comparison. Might cause error if there are missing values.')
    feature <- mmo$feature_data
  } else {
    feature <- mmo$imputed_feature_data
  }
  metadata <- mmo$metadata
  #Get sample names
  group1_samples <- metadata |> dplyr::filter(.data$group == group1) |> dplyr::pull(sample)
  group2_samples <- metadata |> dplyr::filter(.data$group == group2) |> dplyr::pull(sample)
  # Extract numeric matrices (features x samples) — avoids 82K-iteration loop
  g1_mat <- as.matrix(feature[, group1_samples, drop = FALSE])
  g2_mat <- as.matrix(feature[, group2_samples, drop = FALSE])

  m1 <- rowMeans(g1_mat, na.rm = TRUE)
  m2 <- rowMeans(g2_mat, na.rm = TRUE)
  n1 <- rowSums(!is.na(g1_mat))
  n2 <- rowSums(!is.na(g2_mat))

  # Per-row variance without a loop: sweep subtracts row means from each row
  row_var <- function(mat, n) {
    mc <- sweep(mat, 1L, rowMeans(mat, na.rm = TRUE), `-`)
    rowSums(mc^2L, na.rm = TRUE) / pmax(n - 1L, 1L)
  }
  v1 <- row_var(g1_mat, n1)
  v2 <- row_var(g2_mat, n2)

  # Vectorised Welch t-test
  se     <- sqrt(v1 / n1 + v2 / n2)
  t_stat <- (m1 - m2) / se
  # Welch-Satterthwaite degrees of freedom
  df     <- (v1/n1 + v2/n2)^2 /
    ((v1/n1)^2 / pmax(n1 - 1L, 1L) + (v2/n2)^2 / pmax(n2 - 1L, 1L))
  pval   <- 2 * pt(-abs(t_stat), df = pmax(df, 1))
  pval[!is.finite(t_stat)] <- 1   # constant features get p = 1

  log2FC <- log2(m2 / m1)
  padj   <- p.adjust(pval, method = correction)
  #Store in results
  results <- data.frame(
    log2FC= log2FC,
    padj = padj,
    pval = pval
  )
  names(results) <- c(paste(group1, "vs", group2, "log2FC", sep = "_"), paste(group1, "vs", group2, "padj", sep = "_"), paste(group1, "vs", group2, "pval", sep = "_"))
  #Add pairwise results to the mmo
  mmo$pairwise <- cbind(mmo$pairwise, results)
  print(paste(group2, '/', group1, 'comparison was completed'))
  return(mmo)
}


#' Generates lists of DAMs (Differentially Accumulated Metabolites) for each comparison in the \code{mmo}
#'
#' This function generates lists of upregulated and downregulated DAMs for each pairwise comparison in the \code{mmo}.
#' It uses log2 fold change and adjusted p-value thresholds to determine significance.
#' Make sure to run PairwiseComp() for all desired comparisons before using this function.
#'
#' @param mmo The \code{mmo} with pairwise comparison matrix
#' @param fc_cutoff The threshold of log2 fold change to be considered significant (default: 0.5849625, which is log2(1.5))
#' @param pval_cutoff The threshold of adjusted p-value to be considered significant (default: 0.05)
#' @param use_padj Boolean value indicating whether to use adjusted p-value (default: TRUE)
#' @return A list containing two lists: DAMs_up and DAMs_down
#' @export
#' @examplesIf FALSE
#' dams <- GetDAMs(mmo, fc_cutoff = 0.5849625, pval_cutoff = 0.05)
#' dams_up <- dams$DAMs_up
#' dams_down <- dams$DAMs_down
GetDAMs <- function(mmo, fc_cutoff = 0.5849625, pval_cutoff = 0.05, use_padj = TRUE) {
  # Generate the list of comparisons automatically by looking up mmo$pairwise
  comparison_columns <- colnames(mmo$pairwise)
  log2FC_columns <- grep("log2FC", comparison_columns, value = TRUE)
  comparisons <- unique(sub("log2FC", "", log2FC_columns))
  comparisons <- sub("_$", "", comparisons)  # Remove trailing underscore from comparisons
  # Make list of DAMs for up and downregulation for each comparison
  DAMs_up <- list()
  DAMs_down <- list()
  for (comp in comparisons) {
    group1 <- strsplit(comp, "_vs_")[[1]][1]
    group2 <- strsplit(comp, "_vs_")[[1]][2]
    if (use_padj) {
      DAMs_up[[paste(comp, "up", sep = ".")]] <- filter(mmo$pairwise, get(paste(comp, "log2FC", sep = "_")) > fc_cutoff & get(paste(comp, "padj", sep = "_")) < pval_cutoff)$id
      DAMs_down[[paste(comp, "down", sep = ".")]] <- filter(mmo$pairwise, get(paste(comp, "log2FC", sep = "_")) < -fc_cutoff & get(paste(comp, "padj", sep = "_")) < pval_cutoff)$id
    } else {
      DAMs_up[[paste(comp, "up", sep = ".")]] <- filter(mmo$pairwise, get(paste(comp, "log2FC", sep = "_")) > fc_cutoff & get(paste(comp, "pval", sep = "_")) < pval_cutoff)$id
      DAMs_down[[paste(comp, "down", sep = ".")]] <- filter(mmo$pairwise, get(paste(comp, "log2FC", sep = "_")) < -fc_cutoff & get(paste(comp, "pval", sep = "_")) < pval_cutoff)$id
    }
  }
  names(DAMs_up) <- paste(comparisons, "up", sep = ".")
  names(DAMs_down) <- paste(comparisons, "down", sep = ".")
  return(list(DAMs_up = DAMs_up, DAMs_down = DAMs_down))
}



#' Volcano plot for visualizing differential metabolite analysis results
#'
#' This function generates a volcano plot using data from \code{mmo$pairwise} (\code{PairwiseComp(mmo, 'group1', 'group2')} should be specified),
#' highlighting upregulated and downregulated features based on log2 fold change and adjusted p-value
#'
#' @param mmo The \code{mmo} with pairwise comparison matrix
#' @param comp The comparison to visualize, e.g., 'group1_vs_group2
#' @param topk The number of top features to label in the plot (default: 10)
#' @param log2FC_thr The threshold of log2 fold change to be considered significant (default: 1)
#' @param pthr The threshold of adjusted p-value to be considered significant (default: 0.05)
#' @param outdir The output file path for the volcano plot (default: 'volcano.png')
#' @param height The height of the output plot in inches (default: 5)
#' @param width The width of the output plot in inches (default: 5)
#' @param save_output A logical value indicating whether to save the output plot (default: TRUE)
#' @param use_padj A logical value indicating whether to use adjusted p-value (default: TRUE)
#' @return A list containing the volcano plot and the data used to generate it
#' @export
#' @examplesIf FALSE
#' VolcanoPlot(
#'  mmo, comp = 'Control_vs_Treatment1',
#'  topk = 10, log2FC_thr = 1, pthr = 0.05,
#'  outdir = 'volcano_con_tre1.png', height = 5, width = 5
#' )
VolcanoPlot <- function(mmo, comp, topk = 10, log2FC_thr = 1,pthr = 0.05, outdir = 'volcano.png', height = 5, width = 5, save_output = TRUE, use_padj = TRUE){
  .require_pkg("ggrepel")
  VolData <- mmo$pairwise |> dplyr::select(.data$id, all_of(c(paste(comp, 'log2FC', sep = '_'), paste(comp, 'padj', sep = '_'), paste(comp, 'pval', sep = '_'))))
  colnames(VolData) <- c('id', 'log2FC', 'padj', 'pval')
  if(use_padj){ ## BE: revised part (for select whether padj is used)
    VolData <- VolData |>
      mutate(Expression = dplyr::case_when(log2FC >= log2FC_thr & padj <= pthr ~ "Up-regulated", ## BE: revised (to use log2FC_thr and pthr)
                              log2FC <= -log2FC_thr & padj <= pthr ~ "Down-regulated",
                              TRUE ~ "Not significant"))
  } else{ ## BE: when using raw pval...
    VolData <- VolData |>
      mutate(Expression = dplyr::case_when(log2FC >= log2FC_thr & pval <= pthr ~ "Up-regulated", ## BE: revised (to use log2FC_thr and pthr)
                              log2FC <= -log2FC_thr & pval <= pthr ~ "Down-regulated",
                              TRUE ~ "Not significant"))
  }
  VolData$Expression <- factor(VolData$Expression, levels = c("Down-regulated", "Not significant", "Up-regulated"))
  if(use_padj){ ## BE: revised part
    top_features <- dplyr::bind_rows(
      VolData |>
        filter(.data$Expression =='Up-regulated')  |>
        arrange(dplyr::desc(abs(.data$log2FC)), .data$padj) |>
        head(topk),
      VolData |>
        filter(.data$Expression == 'Down-regulated') |>
        arrange(dplyr::desc(abs(.data$log2FC)), .data$padj) |>
        head(topk)
    )
  } else{ ##BE: revised part
    top_features <- dplyr::bind_rows(
      VolData |>
        filter(.data$Expression =='Up-regulated')  |>
        arrange(dplyr::desc(abs(.data$log2FC)), .data$pval) |>
        head(topk),
      VolData |>
        filter(.data$Expression == 'Down-regulated') |>
        arrange(dplyr::desc(abs(.data$log2FC)), .data$pval) |>
        head(topk)
    )
  }

  if(use_padj){ ## BE: revised part
    ## padj ver
    plot <- ggplot(VolData, aes(x = .data$log2FC, y = -log(.data$padj, 10))) +
      geom_point(aes(color = .data$Expression), size = 0.4)+
      xlab(expression("log"[2]*"FC")) + ylab(expression("-log"[10]*"FDR"))+
      scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
      guides(colour = guide_legend(override.aes = list(size=1.5))) +
      theme_classic()+
      ggrepel::geom_label_repel(data = top_features,
                      mapping = aes(.data$log2FC, -log(.data$padj,10), label = .data$id), #previous: label = .data$id
                      size = 2)
  } else{ ## using raw pval
      plot <- ggplot(VolData, aes(x = .data$log2FC, y = -log(.data$pval, 10))) +
      geom_point(aes(color = .data$Expression), size = 0.4)+
      xlab(expression("log"[2]*"FC")) + ylab(expression("-log"[10]*"p-value"))+
      scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
      guides(colour = guide_legend(override.aes = list(size=1.5))) +
      theme_classic()+
      ggrepel::geom_label_repel(data = top_features,
                      mapping = aes(.data$log2FC, -log(.data$padj,10), label = .data$id), #previous: label = .data$id
                      size = 2)
  }

  ## BE edited...
  plot <- plot + ggtitle(paste(comp, "; Up=", sum(VolData$Expression == "Up-regulated"),
                                        "; Down=", sum(VolData$Expression == "Down-regulated"),
                                        "; NS=", sum(VolData$Expression == "Not significant"), sep="")) +
               theme(plot.title = element_text(hjust = 0.0, size=6))

  #plot
  print(plot) ## BE: revised in case of the function call in the for-loop

  if (save_output){
    ggsave(outdir, height = height, width = width)
    readr::write_csv(VolData, paste0(outdir, '_volcano_data.csv'))
  }
  return(list(plot = plot, df = VolData))
}

########################################################################################
# Define functions for multivariate analysis
########################################################################################


#' Plots PCA and performs PERMANOVA
#'
#' This function performs PCA analysis and generates a PCA plot with optional filtering of features and groups.
#' It also conducts PERMANOVA and saves the results to CSV files.
#'
#' @param mmo The \code{mmo} with feature data and metadata
#' @param color A vector of colors for the groups in the plot. Make sure the names correspond to the group names in metadata
#' @param outdir The output file path for the PCA plot (default: 'PCA')
#' @param normalization The normalization method to use for feature data. Options are 'None', 'Log', 'Meancentered', or 'Z' (default: 'Z')
#' @param filter_id Boolean to filter features by id_list (default: FALSE)
#' @param id_list A vector of feature names to filter (default: NULL)
#' @param filter_group Boolean to filter groups by group_list (default: FALSE)
#' @param group_list A vector of group names to filter (default: NULL)
#' @param label Boolean to indicate whether to label points with sample names (default: TRUE)
#' @param save_output Boolean; if TRUE (default) write plot (.pdf) and PERMANOVA
#'   tables using `outdir` as prefix. If FALSE, nothing is written.
#' @return A list with elements `plot` (ggplot), `df` (raw data to generate plots),
#'   and `permanova` (results from `permanova_stat`).
#' @export
#' @examplesIf FALSE
#' PCAplot(
#'  mmo, color = c("Control" = "blue", "Treatment1" = "red", "Treatment2" = "green"),
#'  outdir = 'PCA_plot', normalization = 'None',
#'  filter_id = FALSE, filter_group = FALSE, label = FALSE
#' )
#' PCAplot(
#'  mmo, color = c("Control" = "blue", "Treatment1" = "red"),
#'  outdir = 'PCA_plot', normalization = 'Z',
#'  filter_id = TRUE, id_list = Glucosinolates,
#'  filter_group = TRUE, group_list = c("Control", "Treatment1"), label = TRUE
#' )
PCAplot <- function(mmo, color = NULL, outdir = 'PCA', normalization = 'Z', filter_id = FALSE, id_list = NULL, filter_group = FALSE, group_list = NULL, label = TRUE, save_output = TRUE){
  .require_pkg("ggrepel")
  .require_pkg("stats")
  if (filter_id||filter_group){
    mmo <- filter_mmo(mmo, id_list = id_list, group_list = group_list)
  }
  metadata <- mmo$metadata
  feature <- GetNormFeature(mmo, normalization)
  # If color is NULL, set default colors
  if (is.null(color)) {
    groups <- unique(metadata$group)
    colors <- setNames(grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(groups)), groups)
  }
  # Perform PCA on normalized feature data
  feature_data_pca <- feature[, -(1:2)]
  feature_data_pca <- t(feature_data_pca) # samples as rows, features as columns
  zero_var_cols <- apply(feature_data_pca, 2, var) == 0 # Remove zero columns
  feature_data_pca <- feature_data_pca[, !zero_var_cols]
  pca_res <- stats::prcomp(feature_data_pca, scale. = TRUE)
  pca_df <- as.data.frame(pca_res$x)
  pca_df$group <- metadata$group[match(rownames(pca_df), metadata$sample)]
  # Get PC1 and PC2 loadings
  pca_res_summary <- summary(pca_res)
  PC1.loading <- as.data.frame(pca_res_summary$importance)$PC1[2]
  PC2.loading <- as.data.frame(pca_res_summary$importance)$PC2[2]

  if (label == TRUE){
    plot <- ggplot(pca_df, aes(x = .data$PC1, y = .data$PC2, color = .data$group, label = rownames(pca_df))) +
      geom_point(size = 3) +
      ggrepel::geom_label_repel(aes(label = rownames(pca_df)), size = 3) +
      theme_classic() +
      labs(x = paste("PC1 (", PC1.loading ,")", sep=""),
           y = paste("PC2 (", PC2.loading ,")", sep="")) +
      scale_color_manual(values = color)+
      stat_ellipse(aes(group = .data$group), level = 0.90)
  } else {
    plot <- ggplot(pca_df, aes(x = .data$PC1, y = .data$PC2, color = .data$group)) +
      geom_point(size = 3) +
      # geom_label_repel(aes(label = rownames(pca_df)), size = 3) +
      theme_classic() +
      labs(x = paste("PC1 (", PC1.loading ,")", sep=""),
           y = paste("PC2 (", PC2.loading ,")", sep="")) +
      scale_color_manual(values = color)+
      stat_ellipse(aes(group = .data$group), level = 0.90)
  }
  permanova <- permanova_stat(feature_data_pca, metadata, mode = 'data', filter_group = filter_group, group_list = group_list)
  plot
  if (isTRUE(save_output)) {
    ggsave(paste0(outdir, '.pdf'), width = 6, height = 6)
    readr::write_csv(permanova$permanova_res, paste0(outdir, '_permanova_results.csv'))
    readr::write_csv(as.data.frame(permanova$pairwise_raw), paste0(outdir, '_pairwise_permanova_results.csv'))
    readr::write_csv(as.data.frame(permanova$pairwise_p_matrix), paste0(outdir, '_pairwise_permanova_pvalue_matrix.csv'))
    readr::write_csv(as.data.frame(permanova$pairwise_F_matrix), paste0(outdir, '_pairwise_permanova_Fvalue_matrix.csv'))
    readr::write_csv(as.data.frame(permanova$pairwise_R2_matrix), paste0(outdir, '_pairwise_permanova_R2_matrix.csv'))
  }

  return(invisible(list(plot = plot, df = pca_df, permanova = permanova)))
}




#' PLS-DA plot with feature loadings
#'
# This function performs PLS-DA analysis and generates a PLS-DA plot with feature loadings.
#'
#' @param mmo The \code{mmo} with feature data and metadata
#' @param color A vector of colors for the groups in the plot. Make sure the names correspond to the group names in metadata
#' @param topk The number of top features to display in the plot (default: 10)
#' @param outdir The output file path for the PLS-DA plot
#' @param normalization The normalization method to use for feature data. Options are 'None', 'Log', 'Meancentered', or 'Z' (default: 'Z')
#' @param filter_id Boolean to filter features by id_list (default: FALSE)
#' @param id_list A vector of feature names to filter (default: NULL)
#' @param filter_group Boolean to filter groups by group_list (default: FALSE)
#' @param group_list A vector of group names to filter (default: NULL)
#' @param save_output Boolean; if TRUE (default) write plot (.pdf) and loadings
#'   tables using `outdir` as prefix. If FALSE, nothing is written.
#' @return A list with elements `plot` (ggplot), `df` (raw data to generate plots),
#'   and `loadings` (loadings for PLSDA).
#' @export
#' @examplesIf FALSE
#' PLSDAplot(
#'  mmo, color = c("Control" = "blue", "Treatment1" = "red", "Treatment2" = "green"),
#'  topk = 10, outdir = 'PLSDA_plot.pdf', normalization = 'Z',
#'  filter_id = FALSE, filter_group = FALSE
#' )
#' PLSDAplot(
#'  mmo, color = c("Control" = "blue", "Treatment1" = "red"),
#'  topk = 5, outdir = 'PLSDA_plot.pdf', normalization = 'Log',
#'  filter_id = TRUE, id_list = Glucosinolates,
#'  filter_group = TRUE, group_list = c("Control", "Treatment1")
#' )
PLSDAplot <- function(mmo, color, topk = 10, outdir, normalization = 'Z', filter_id = FALSE, id_list = NULL, filter_group = FALSE, group_list = NULL, save_output = TRUE) {
  .require_pkg("caret")
  .require_pkg("ggrepel")
  if (filter_id||filter_group){
    mmo <- filter_mmo(mmo, id_list = id_list, group_list = group_list)
  }
  metadata <- mmo$metadata
  #Get appropriate feature by normalization parameter
  feature <- GetNormFeature(mmo, normalization)

  # All feature or filtered feature
  X <- t(as.matrix(feature[, -(1:2)]))
  Y <- c()
  for (col in colnames(feature)[-c(1,2)]){
    Y <- append(Y, metadata[metadata$sample == col, ]$group)
  }
  Y <- as.factor(Y)

  plsda_model <- caret::plsda(X, Y, ncomp = 2)
  scores <- plsda_model$scores[, 1:2]
  plsda_df <- data.frame(Comp1 = scores[, 1], Comp2 = scores[, 2], Group = Y)
  loadings <- plsda_model$loadings
  loadings_comp1 <- loadings[, 1]
  loadings_comp2 <- loadings[, 2]
  loadings_df <- data.frame(ID = mmo$feature_data$id,
                            Comp1_Loading = loadings_comp1,
                            Comp2_Loading = loadings_comp2)

  top_features <- loadings_df |>
  mutate(abs_loading_comp1 = abs(.data$Comp1_Loading),
         abs_loading_comp2 = abs(.data$Comp2_Loading)) |>
  arrange(dplyr::desc(.data$abs_loading_comp1 + .data$abs_loading_comp2)) |>
  head(topk)
  loading_scale <- 1
  if (topk > 0){
  loading_scale <- max(abs(scores))/(4*max(abs(top_features$Comp1_Loading)))
  }

 plot <-  ggplot(plsda_df, aes(x = .data$Comp1, y = .data$Comp2, color = .data$Group)) +
    geom_point(size = 3) +
    theme_classic() +
    stat_ellipse(level = 0.90) +
    ggtitle("PLS-DA Plot") +
    labs(x = "Component 1", y = "Component 2") +
    scale_color_manual(values = color) +
    theme(legend.position = "right")+
    geom_segment(data = top_features,
                aes(x = 0, y = 0, xend = .data$Comp1_Loading * loading_scale, yend = .data$Comp2_Loading * loading_scale),  # Scale the arrows
                arrow = grid::arrow(length = grid::unit(0.3, "cm")), color = "grey", linewidth = 1) +
    # Add labels for the top 10 features
    ggrepel::geom_text_repel(data = top_features,
              aes(x = .data$Comp1_Loading * loading_scale, y = .data$Comp2_Loading * loading_scale, label = .data$ID),
              color = "black", vjust = 1.5, size = 3)
  plot
  if (save_output) {
    ggsave(paste0(outdir, 'PLSDA_plot.pdf'), height = 6, width = 6)
    readr::write_csv(loadings_df, paste0(outdir, 'PLSDA_loadings.csv'))
  }
  return(list(plot = plot, df = plsda_df, loadings = loadings_df))
}

#' Generate input files to be used for pheatmap from the \code{mmo}
#'
#' This function generates heatmap inputs from the \code{mmo}, including fold change or mean values,
#' distance matrix, and row labels for custom-annotated features.
#'
#' @param mmo The \code{mmo} with sirius annotation and normalized data
#' @param filter_id Boolean to filter features by id_list (default: FALSE)
#' @param id_list A vector of feature names to filter (default: NULL)
#' @param filter_group Boolean to filter groups by group_list (default: FALSE)
#' @param group_list A vector of group names to filter (default: NULL)
#' @param summarize The summarization method to use. Options are 'fold_change' or 'mean' (default: 'mean')
#' @param control_group The group to use as control for fold change calculation (default: 'ctrl')
#' @param normalization The normalization method to use. Options are 'None', 'Log', 'Meancentered', or 'Z'
#' @param distance The distance metric to use. Options are 'dreams', 'cosine', or 'm2ds' (default: 'dreams')
#' @return A list containing the following elements:
#' - FC_matrix: A matrix of fold change or mean values
#' - dist_matrix: A distance matrix based on the specified distance metric
#' - row_label: A vector of row labels for custom-annotated features (See AddCustomAnnot()). If no custom annotation is available, feature IDs are used.
#' - heatmap_data: A data frame containing the heatmap data with feature IDs and values
#' @export
#' @examplesIf FALSE
#' # Generate heatmap inputs to visualize fold change values with log normalization and dreams distance
#' heatmap_inputs <- GenerateHeatmapInputs(
#'  mmo, summarize = 'fold_change', control_group = 'Control',
#'  normalization = 'None', distance = 'dreams'
#' )
#' # Generate heatmap inputs to visualize mean values
#' heatmap_inputs <- GenerateHeatmapInputs(
#'  mmo, summarize = 'mean', normalization = 'None', distance = 'dreams'
#' )
#' # The resulting list contains FC_matrix, dist_matrix, row_label, and heatmap_data
#' # A heatmap can be generated using pheatmap
#' # 'clustering_distance_rows' option make the dendrogram follows chemical distances of features.
#' #  -Delete this option to visualize the heatmap following cannonical clustering
#' pheatmap(mat = heatmap_inputs$FC_matrix,
#'     cluster_rows = TRUE, #do not change
#'     clustering_distance_rows = heatmap_inputs$dist_matrix,
#'     cluster_cols = TRUE,
#'     clustering_method = "average", #UPGMA
#'     show_rownames = TRUE,
#'     show_colnames = TRUE,
#'     cellwidth = 25,
#'     cellheight = 0.05,
#'     treeheight_row = 100,
#'     fontsize_row = 3,
#'     fontsize_col = 15,
#'     scale = 'none',
#'     annotation_names_row = TRUE,
#'     labels_row = heatmap_inputs$row_label,
#'     )
# DEPRECATED: superseded by GenerateHeatmapInputs().
# GenerateHeatmapInputs_derep() uses GetDistanceMat() (.dissim slots).
# GenerateHeatmapInputs() uses GetSimMat() (.sim slots) and converts
# similarity to distance (1-S) when building the dist object for clustering.
GenerateHeatmapInputs_derep <- function(mmo, filter_id = FALSE, id_list = NULL,
                                filter_group = FALSE, group_list = NULL,
                                summarize = 'mean', control_group = 'ctrl',
                                normalization = 'None', distance = NULL) {
  if (filter_id||filter_group){
    mmo <- filter_mmo(mmo, id_list = id_list, group_list = group_list)
  }
  # 12.1.1. Get summarized data (group mean or FC)
  group_means <- GetGroupMeans(mmo, normalization = normalization)
  if (summarize == 'fold_change'){
    fold_change <- GetLog2FoldChange(group_means, control_group = control_group)
    heatmap_data <- fold_change
    heatmap_data[[control_group]] <- NULL
  } else if(summarize == 'mean'){
    heatmap_data <- group_means
  }
  # 12.1.2. Filter features
  # Determine distance metric
  if (!is.null(distance)){
    distance_matrix <- GetDistanceMat(mmo, distance = distance)
    heatmap_data <- heatmap_data |> filter(.data$id %in% rownames(distance_matrix)) # remove features not in distance matrix
  }

  # make matrix for heatmap
  FC_matrix <- as.matrix(heatmap_data[,-1])
  rownames(FC_matrix) <- heatmap_data$id
  if (!is.null(distance)){
    # Reorder the rows of distance_matrix to match the order of FC_matrix_
    distance_matrix <- distance_matrix[rownames(FC_matrix), rownames(FC_matrix)]
    dist_matrix <- as.dist(distance_matrix)
  } else{
    dist_matrix <- NULL
  }

  row_label <- rownames(FC_matrix)
  return(list(FC_matrix = FC_matrix, dist_matrix = dist_matrix, row_label = row_label, heatmap_data = heatmap_data))
}


#' GenerateHeatmapInputs (similarity-based)
#'
#' Identical to \code{GenerateHeatmapInputs_derep()} but retrieves feature
#' similarity matrices via \code{GetSimMat()} (.sim slots) and converts to
#' distance (\code{1 - S}) when building the \code{dist} object for heatmap
#' row clustering.
#'
#' @inheritParams GenerateHeatmapInputs_derep
#' @export
GenerateHeatmapInputs <- function(mmo, filter_id = FALSE, id_list = NULL,
                                  filter_group = FALSE, group_list = NULL,
                                  summarize = 'mean', control_group = 'ctrl',
                                  normalization = 'None', distance = NULL) {
  if (filter_id||filter_group){
    mmo <- filter_mmo(mmo, id_list = id_list, group_list = group_list)
  }
  # 12.1.1. Get summarized data (group mean or FC)
  group_means <- GetGroupMeans(mmo, normalization = normalization)
  if (summarize == 'fold_change'){
    fold_change <- GetLog2FoldChange(group_means, control_group = control_group)
    heatmap_data <- fold_change
    heatmap_data[[control_group]] <- NULL
  } else if(summarize == 'mean'){
    heatmap_data <- group_means
  }
  # 12.1.2. Filter features; use similarity matrix from .sim slot
  if (!is.null(distance)){
    sim_matrix <- GetSimMat(mmo, distance = distance)
    heatmap_data <- heatmap_data |> filter(.data$id %in% rownames(sim_matrix))
  }

  # make matrix for heatmap
  FC_matrix <- as.matrix(heatmap_data[,-1])
  rownames(FC_matrix) <- heatmap_data$id
  if (!is.null(distance)){
    # Reorder similarity matrix to match FC_matrix row order
    sim_matrix <- sim_matrix[rownames(FC_matrix), rownames(FC_matrix)]
    # A heatmap needs a dense distance matrix, so this cannot stay sparse.
    # Stop before allocating tens of gigabytes rather than exhausting memory.
    n_keep <- nrow(sim_matrix)
    if (n_keep > 10000L) {
      stop("GenerateHeatmapInputs: clustering ", n_keep, " features requires a dense ",
           n_keep, " x ", n_keep, " distance matrix (about ",
           format(round(n_keep^2 * 8 / 2^30, 1), nsmall = 1), " GB).\n",
           "Use filter_id = TRUE with id_list to select fewer features, or pass ",
           "distance = NULL to skip chemical clustering.", call. = FALSE)
    }
    # Convert similarity to distance (1 - S) for heatmap clustering
    dist_matrix <- as.dist(1 - as.matrix(sim_matrix))
  } else{
    dist_matrix <- NULL
  }

  row_label <- rownames(FC_matrix)
  return(list(FC_matrix = FC_matrix, dist_matrix = dist_matrix, row_label = row_label, heatmap_data = heatmap_data))
}

#' PlotNPCStackedBar
#'
#' This function generates a stacked bar plot showing the count of features in each group categorized by NPC_pathway.
#' It uses the \code{mmo} with sirius annotation and normalized data.
#' Make sure you don't run ReplaceZero() before using this function, as it may remove presence/absence information.
#'
#' @param mmo The \code{mmo} with sirius annotation and normalized data
#' @param group_col The column name in metadata to use for grouping samples
#' @param outdir The output file path for the stacked bar plot (e.g., 'NPC_stacked_bar.png')
#' @param width The width of the output plot
#' @param height The height of the output plot
#' @param save_output boolean, whether to save the output plot
#' @return A list containing the stacked bar plot and the data used to generate it
#' @export
#' @examplesIf FALSE
#' PlotNPCStackedBar(
#'  mmo, group_col = 'treatment',
#'  outdir = 'NPC_stacked_bar.png', width = 6, height = 3
#' )
PlotNPCStackedBar <- function(mmo, group_col, outdir, width = 6, height = 3, save_output = TRUE) {
  mmo <- SwitchGroup(mmo, group_col)
  feature_data <- mmo$feature_data
  metadata <- mmo$metadata

  # For each group, get features present in any sample
  group_features <- lapply(unique(metadata$group), function(grp) {
    samples <- metadata |> filter(.data$group == grp) |> pull(.data$sample)
    present <- feature_data |> dplyr::select(all_of(samples))
    feature_data$feature[base::rowSums(!is.na(present) & present > 0) > 0]
  })
  names(group_features) <- unique(metadata$group)

  # Build long data frame: feature, group, NPC_pathway
  annot <- mmo$sirius_annot[, c("feature", "NPC#pathway")]
  colnames(annot) <- c("feature", "NPC_pathway")

  bar_df <- do.call(rbind, lapply(names(group_features), function(grp) {
    data.frame(
      feature = group_features[[grp]],
      group = grp,
      stringsAsFactors = FALSE
    )
  }))
  bar_df <- merge(bar_df, annot, by = "feature", all.x = TRUE)
  bar_df <- bar_df[!is.na(bar_df$NPC_pathway) & bar_df$NPC_pathway != "", ]

  # Count features per group and NPC_pathway
  plot_df <- bar_df |>
    dplyr::group_by(.data$group, .data$NPC_pathway) |>
    dplyr::summarise(count = dplyr::n(), .groups = "drop")

  # Set colors for NPC_pathway
  npc_pathways <- unique(plot_df$NPC_pathway)
  .require_pkg("RColorBrewer")
  bar_colors <- setNames(RColorBrewer::brewer.pal(min(length(npc_pathways), 8), "Set2"), npc_pathways)

  # Plot stacked bar by NPC_pathway
  plot <- ggplot(plot_df, aes(x = .data$group, y = .data$count, fill = .data$NPC_pathway, label = .data$count)) +
    geom_bar(stat = "identity", position = "stack") +
    geom_text(aes(group = .data$NPC_pathway), position = position_stack(vjust = 0.5), size = 3, color = "white", fontface = "bold") +
    scale_fill_manual(values = bar_colors) +
    coord_flip() +
    labs(x = "Group", y = "Feature Count", fill = "NPC Pathway") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  if (save_output){
    ggsave(paste0(outdir, 'NPC_stacked_bar.pdf'), width = width, height = height)
    readr::write_csv(plot_df, paste0(outdir, 'NPC_stacked_bar.csv'))
  }
  return(list(plot = plot, df = plot_df))
}

#' Enrichment analysis for Canopus-predicted terms
#'
#' This function performs enrichment analysis for Canopus-predicted terms on a given list of features.
#'
#' @param mmo The \code{mmo} with sirius annotation and normalized data
#' @param list_test A vector containing ids of features to analyze
#' @param pthr The threshold for adjusted p-value to be considered significant (default: 0.1)
#' @param sig A logical value indicating whether to return only significant terms (default: TRUE)
#' @param term_level The level of term to use for enrichment analysis
#'               Options are 'NPC_pathway', 'NPC_superclass', 'NPC_class', 'ClassyFire_superclass', 'ClassyFire_class',
#'              'ClassyFire_subclass', 'ClassyFire_level5', or 'ClassyFire_most_specific' (default: 'NPC_pathway')
#' @param representation The representation type for enrichment analysis. Options are 'greater' for overrepresentation (default: 'greater')
#' @param pval pvalue options-pval or fdr (default: 'pval')
#' @return A data frame containing the enrichment results, including term level,
#'        term name, subset count, total count, fold enrichment, p-value, and adjusted p-value (FDR)
#' @export
#' @examplesIf FALSE
#' # Perform enrichment analysis for a list of features using NPC_pathway level
#' sig_terms <- CanopusLevelEnrichmentAnal(
#'  mmo, list_test = c("feature1", "feature2"), pthr = 0.1,
#'  sig = TRUE, term_level = 'NPC_pathway', representation = 'greater'
#' )
#' # Perform enrichment analysis for a list of features using
#' # ClassyFire_class level and return all terms
#' all_terms <- CanopusLevelEnrichmentAnal(
#'  mmo, list_test = c("feature1", "feature2"), pthr = 0.1,
#'  sig = FALSE, term_level = 'ClassyFire_class', representation = 'greater'
#' )
CanopusLevelEnrichmentAnal <- function(mmo,list_test, pthr = 0.1, sig=TRUE, term_level = 'NPC_pathway', representation = 'greater', pval = 'pval'){
  all_feature <- mmo$sirius_annot
  subset_feature <- mmo$sirius_annot |> filter(.data$id %in% list_test)
  # print(paste('total features:', nrow(all_feature), 'list_test features:', nrow(subset_feature)))
  # Select the appropriate term level for enrichment analysis
  if (term_level == "NPC_pathway") {
    all_feature$classifications_split <- all_feature[['NPC#pathway']]
    subset_feature$classifications_split <- subset_feature[['NPC#pathway']]
  } else if (term_level == "NPC_superclass") {
    all_feature$classifications_split <- all_feature[['NPC#superclass']]
    subset_feature$classifications_split <- subset_feature[['NPC#superclass']]
  } else if (term_level == "NPC_class") {
    all_feature$classifications_split <- all_feature[['NPC#class']]
    subset_feature$classifications_split <- subset_feature[['NPC#class']]
  } else if (term_level == "ClassyFire_superclass") {
    all_feature$classifications_split <- all_feature[['ClassyFire#superclass']]
    subset_feature$classifications_split <- subset_feature[['ClassyFire#superclass']]
  } else if (term_level == "ClassyFire_class") {
    all_feature$classifications_split <- all_feature[['ClassyFire#class']]
    subset_feature$classifications_split <- subset_feature[['ClassyFire#class']]
  } else if (term_level == "ClassyFire_subclass") {
    all_feature$classifications_split <- all_feature[['ClassyFire#subclass']]
    subset_feature$classifications_split <- subset_feature[['ClassyFire#subclass']]
  } else if (term_level == "ClassyFire_level5") {
    all_feature$classifications_split <- all_feature[['ClassyFire#level 5']]
    subset_feature$classifications_split <- subset_feature[['ClassyFire#level 5']]
  } else if (term_level == "ClassyFire_most_specific") {
    all_feature$classifications_split <- all_feature[['ClassyFire#most specific class']]
    subset_feature$classifications_split <- subset_feature[['ClassyFire#most specific class']]
  } else {
    stop("Invalid term level. Please choose a valid term level.")
  }

  total_term_counts <- table(unlist(all_feature$classifications_split))
  subset_term_counts <- table(unlist(subset_feature$classifications_split))

  total_term_counts['None'] <- sum(is.na(all_feature$classifications_split))
  subset_term_counts['None'] <- sum(is.na(subset_feature$classifications_split))

  # Perform enrichment analysis using Fisher's exact test
  enrichment_results <- sapply(names(subset_term_counts), function(term) {
    contingency_matrix <- matrix(c(
      subset_term_counts[[term]],
      sum(subset_term_counts) - subset_term_counts[[term]],
      total_term_counts[[term]],
      sum(total_term_counts) - total_term_counts[[term]]
    ), nrow = 2, byrow = TRUE)
    fisher.test(contingency_matrix, alternative = representation)$p.value
  })
  # Adjust p-values for multiple testing
  adjusted_pvalues <- p.adjust(enrichment_results, method = "fdr")
  # Create a results dataframe
  results <- data.frame(
    term_level = term_level,
    term = names(enrichment_results),
    subsetcount = as.numeric(subset_term_counts[names(enrichment_results)]),
    totalcount = as.numeric(total_term_counts[names(enrichment_results)]),
    foldenrichment = (as.numeric(subset_term_counts[names(enrichment_results)]) / nrow(subset_feature))/(as.numeric(total_term_counts[names(enrichment_results)]) / nrow(all_feature)),
    pval = enrichment_results,
    fdr = adjusted_pvalues
  )
  results <- results |> filter(.data$term != 'None')
  # Filter for significantly enriched terms
  if(pval == 'pval'){
    significant_terms <- results |>
      filter(.data$pval < pthr) |>
      arrange(.data$pval)
  } else if (pval == 'fdr'){
    significant_terms <- results |>
      filter(.data$fdr < pthr) |>
      arrange(.data$fdr)
  } else {
    stop("Invalid pval option. Please choose 'pval' or 'fdr'.")
  }
  if(sig==TRUE){
    return(significant_terms)
  }else{
    return(results)
  }
}

#' Generate a plot for enrichment analysis of Canopus-predicted terms
#'
#' This function generates a plot for enrichment analysis of Canopus-predicted terms,
#' showing fold enrichment, p-value, and subset count for each term level.
#'
#' @param mmo The \code{mmo} with sirius annotation and normalized data
#' @param id_list A vector containing names of features to analyze
#' @param pthr The threshold for adjusted p-value to be considered significant (default: 0.05)
#' @param outdir The output file path for the enrichment plot
#' @param height The height of the output plot in inches (default: 5)
#' @param width The width of the output plot in inches (default: 5)
#' @param pval pvalue options-pval or fdr (default: 'pval')
#' @param save_output boolean, whether to save the output plot (default: TRUE)
#' @return A list containing the enrichment plot and the enrichment results
#' @export
#' @examplesIf FALSE
#' CanopusListEnrichmentPlot(
#'  mmo, id_list = DAMs_up$control_vs_treatment1.up,
#'  pthr = 0.05, outdir = 'canopus_enrichment_plot.pdf',
#'  height = 5, width = 5
#' )
#'
CanopusListEnrichmentPlot <- function(mmo, id_list, pthr = 0.05, outdir, height = 5, width = 5, pval = 'pval', save_output = TRUE){
  term_levels = c('NPC_class', 'NPC_superclass', 'NPC_pathway', 'ClassyFire_superclass', 'ClassyFire_class', 'ClassyFire_subclass', 'ClassyFire_level5', 'ClassyFire_most_specific')
  sig.canopus <- data.frame(term = character(),  term_level = character(),subsetcount = double(), totalcount = double(), foldenrichment = double(), pval = double(), fdr = double())
  for (term_level in term_levels){
    sig.canopus <- rbind(sig.canopus, CanopusLevelEnrichmentAnal(mmo, id_list, pthr = pthr, sig = TRUE, term_level = term_level, representation = 'greater', pval = pval))
  }
  sig.canopus <- sig.canopus |> arrange(dplyr::desc(.data$foldenrichment))
  plot <- ggplot(sig.canopus, aes(x = .data$foldenrichment, y = reorder(.data$term, .data$foldenrichment), color = -log(.data$fdr), size = .data$subsetcount)) +
    geom_point() +
    scale_color_gradient(low = 'grey', high = 'red') +
    theme_classic()+
    facet_grid(term_level ~ ., scales = 'free_y', space = 'free', switch = 'y')+
    xlim(0,max(sig.canopus$foldenrichment+1))
    #facet_wrap(~term_level, ncol = 1, scales = 'free_y', strip.position = 'right', shrink = TRUE)
  plot

  if (save_output){
    ggsave(paste0(outdir, 'enrichment.pdf'), height = height, width = width)
    readr::write_csv(sig.canopus, paste0(outdir, 'enrichment.csv'))
  }
  return(list(plot = plot, df = sig.canopus))
}

#' Generate a plot for enrichment analysis of Canopus-predicted terms across multiple levels
#'
#' This function generates a plot for enrichment analysis of Canopus-predicted terms across multiple levels,
#' showing fold enrichment, p-value, and subset count for each term level.
#'
#' @param mmo The \code{mmo} with sirius annotation and normalized data
#' @param id_list A vector containing names of features to analyze
#' @param pthr The threshold for adjusted p-value to be considered significant (default: 0.05)
#' @param outdir The output file path for the enrichment plot
#' @param height The height of the output plot in inches (default: 5)
#' @param width The width of the output plot in inches (default: 5)
#' @param topn The number of top terms to display in the plot (default: 5)
#' @param pval pvalue options-pval or fdr (default: 'pval')
#' @param save_output boolean, whether to save the output plot (default: TRUE)
#' @return A list containing the enrichment plot and the enrichment results
#' @export
#' @examplesIf FALSE
#' CanopusListEnrichmentPlot_2(
#'  mmo, id_list = DAMs_up$control_vs_treatment1.up,
#'  pthr = 0.05, outdir = 'canopus_enrichment_plot_topn.pdf',
#'  height = 5, width = 5, topn = 5
#' )
CanopusListEnrichmentPlot_2 <- function(mmo, id_list, pthr = 0.05, outdir, height = 5, width = 5, topn = 5, pval = 'pval', save_output = TRUE){
  term_levels = c('NPC_class', 'NPC_superclass', 'NPC_pathway', 'ClassyFire_superclass', 'ClassyFire_class', 'ClassyFire_subclass', 'ClassyFire_level5', 'ClassyFire_most_specific')
  sig.canopus <- data.frame(term = character(),  term_level = character(),subsetcount = double(), totalcount = double(), foldenrichment = double(), pval = double(), fdr = double())
  for (term_level in term_levels){
    sig.canopus <- rbind(sig.canopus, CanopusLevelEnrichmentAnal(mmo, id_list, pthr = pthr, sig = TRUE, term_level = term_level, representation = 'greater', pval = pval))
  }
  sig.canopus$term <- paste(sig.canopus$term, ';', sig.canopus$term_level)
  sig.canopus <- sig.canopus |> dplyr::slice_max(order_by = -.data$pval, n = topn)
  sig.canopus <- sig.canopus |> dplyr::arrange(dplyr::desc(.data$foldenrichment))
  plot <- ggplot(sig.canopus, aes(x = .data$foldenrichment, y = reorder(.data$term, .data$foldenrichment), color = -log(.data$fdr), size = .data$subsetcount)) +
    geom_point() +
    scale_color_gradient(low = 'grey', high = 'red') +
    theme_classic()+
    xlim(0,max(sig.canopus$foldenrichment+1))+
    ylab('Chemical Class')
  plot
  if (save_output){
    ggsave(outdir, height = height, width = width)
    readr::write_csv(sig.canopus, paste0(outdir, 'enrichment.csv'))
  }
  return(list(plot = plot, df = sig.canopus))
}

#' Generate a plot for enrichment analysis of Canopus-predicted terms at a specific level using a list of vectors of features
#'
#' This function generates a plot for enrichment analysis of Canopus-predicted terms at a specific level,
#' showing fold enrichment, p-value, and subset count for each term.
#'
#' @param mmo The \code{mmo} with sirius annotation and normalized data
#' @param comp.list A list to analyze, where each element is a vector of feature ids
#' @param term_level The level of term to use for enrichment analysis.
#'               Options are 'NPC_pathway', 'NPC_superclass', 'NPC_class',
#'              'ClassyFire_superclass', 'ClassyFire_class', 'ClassyFire_subclass',
#'             'ClassyFire_level5', or 'ClassyFire_most_specific' (default: 'NPC_pathway')
#' @param pthr The threshold for adjusted p-value to be considered significant (default: 0.1)
#' @param representation The representation type for enrichment analysis. Options are 'greater' for overrepresentation (default: 'greater')
#' @param outdir The output directory for saving the plot and the enrichment results (default: 'enrichment')
#' @param height The height of the output plot in inches (default: 5)
#' @param width The width of the output plot in inches (default: 5)
#' @param pval pvalue options-pval or fdr (default: 'pval')
#' @param save_output boolean, whether to save the output plot (default: TRUE)
#' @return A list containing the enrichment plot and the enrichment results
#' @export
#' @examplesIf FALSE
#' # Perform enrichment analysis for multiple comparisons using NPC_pathway level
#' comp.list <- list(
#'   comparison1 = DAMs_up$control_vs_treatment1.up,
#'   comparison2 = DAMs_up$control_vs_treatment2.up
#' )
#' CanopusLevelEnrichmentPlot(
#'  mmo, comp.list = comp.list, term_level = 'NPC_pathway',
#'  pthr = 0.1, representation = 'greater', outdir = 'enrichment_plot',
#'  height = 5, width = 5
#' )
CanopusLevelEnrichmentPlot <- function(mmo = mmo, comp.list, term_level = 'NPC_pathway',pthr = 0.1, representation = 'greater', outdir = 'enrichment', height = 5, width = 5, pval = 'pval', save_output = TRUE){
  df.EA <- data.frame()
  sig.terms <- c()
  for(list in names(comp.list)){
    # Calculate enrichment score for all terms
    res <- CanopusLevelEnrichmentAnal(mmo = mmo, list_test = comp.list[[list]],sig=FALSE, pthr = pthr, representation = representation, term_level = term_level, pval = pval)
    res <- res |> mutate(comp = list)
    df.EA <- dplyr::bind_rows(df.EA, res)
    # get terms that are at least once enriched in one comparison
    res.sig <- CanopusLevelEnrichmentAnal(mmo = mmo, list_test = comp.list[[list]],sig=TRUE, pthr = pthr, representation = representation, term_level = term_level, pval = pval)
    sig.terms <- append(sig.terms, res.sig$term)
  }
  sig.terms <- unique(sig.terms)
  df.EA.sig <- df.EA |> filter(.data$term %in% sig.terms)
  if (pval == 'pval'){
    df.EA.sig <- df.EA.sig |>
      mutate(label = cut(
        .data$pval,
        breaks = c(0,0.001, 0.01, 0.05, 0.1, 1),
        labels = c("***", "**", "*", ".", "")
      ))
    plot <- ggplot(data = df.EA.sig, aes(x = .data$comp, y = .data$term, label = .data$label))+
      geom_point(aes(size = .data$subsetcount, color = .data$pval))+
      geom_text()+
      scale_size_area(name = 'Count', max_size = 10)+
      scale_color_gradient2(low = 'red', high = 'grey', mid = 'grey', midpoint = 0.4)+
      theme_minimal()+
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))+
      xlab('Comparisons')+
      ylab('Chemical classes')
  } else if (pval == 'fdr'){
    df.EA.sig <- df.EA.sig |>
      mutate(label = cut(
        .data$fdr,
        breaks = c(0,0.001, 0.01, 0.05, 0.1, 1),
        labels = c("***", "**", "*", ".", "")
      ))
    plot <- ggplot(data = df.EA.sig, aes(x = .data$comp, y = .data$term, label = .data$label))+
      geom_point(aes(size = .data$subsetcount, color = .data$fdr))+
      geom_text()+
      scale_size_area(name = 'Count', max_size = 10)+
      scale_color_gradient2(low = 'red', high = 'grey', mid = 'grey', midpoint = 0.4)+
      theme_minimal()+
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))+
      xlab('Comparisons')+
      ylab('Chemical classes')
  }
  plot
  if (save_output){
    ggsave(paste0(outdir, '.pdf'), width = width, height = height)
    readr::write_csv(df.EA, paste0(outdir, '.csv'))
    readr::write_csv(df.EA.sig, paste0(outdir, '_sig.csv'))
  }
  return(list(plot = plot, df = df.EA))
}

#' Generate a plot for enrichment analysis of Canopus-predicted terms across all levels
#'
#' This function generates a plot for enrichment analysis of Canopus-predicted terms across all levels,
#' showing fold enrichment, p-value, and subset count for each term level.
#' @param mmo The \code{mmo} with sirius annotation and normalized data
#' @param comp.list A list to analyze, where each element is a vector of feature names
#' @param terms The terms to analyze. Options are 'all_terms', 'NPC', 'ClassyFire', or 'custom' (default: 'all_terms')
#' @param term_levels list of custom term levels to use
#' @param pthr The threshold for adjusted p-value to be considered significant (default: 0.1)
#' @param representation The representation type for enrichment analysis. Options are 'greater' for overrepresentation (default: 'greater')
#' @param outdir The output directory for saving the plot and the enrichment results (default: 'enrichment')
#' @param height The height of the output plot in inches (default: 10)
#' @param width The width of the output plot in inches (default: 8)
#' @param pval pvalue options-pval or fdr (default: 'pval')
#' @param save_output boolean, whether to save the output plot (default: TRUE)
#' @return A list of the plot and the enrichment results
#' @export
#' @examplesIf FALSE
#' comp.list <- list(
#'   comparison1 = DAMs_up$control_vs_treatment1.up,
#'   comparison2 = DAMs_up$control_vs_treatment2.up
#' )
#' CanopusAllLevelEnrichmentPlot(
#'  mmo, comp.list = comp.list, terms = 'all_terms',
#'  pthr = 0.1, representation = 'greater', outdir = 'enrichment_all_levels',
#'  height = 10, width = 8
#' )
#' CanopusAllLevelEnrichmentPlot(
#'  mmo, comp.list = comp.list, terms = 'NPC',
#'  pthr = 0.1, representation = 'greater', outdir = 'enrichment_NPC_levels',
#'  height = 10, width = 8
#' )
#' CanopusAllLevelEnrichmentPlot(
#'  mmo, comp.list = comp.list, terms = 'ClassyFire',
#'  pthr = 0.1, representation = 'greater', outdir = 'enrichment_ClassyFire_levels',
#'  height = 10, width = 8
#' )
CanopusAllLevelEnrichmentPlot <- function(mmo = mmo, comp.list, terms = 'all_terms', term_levels = NULL, pthr = 0.1, representation = 'greater', outdir = 'enrichment', height = 10, width = 8, pval = 'pval', save_output = TRUE){
  df.EA <- data.frame()
  sig.terms <- c()
  if(terms == 'all_terms'){
    term_levels = c('NPC_class', 'NPC_superclass', 'NPC_pathway', 'ClassyFire_superclass', 'ClassyFire_class', 'ClassyFire_subclass', 'ClassyFire_level5', 'ClassyFire_most_specific')
  } else if (terms == 'NPC'){
    term_levels = c('NPC_class', 'NPC_superclass', 'NPC_pathway')
  } else if (terms == 'ClassyFire'){
    term_levels = c('ClassyFire_superclass', 'ClassyFire_class', 'ClassyFire_subclass', 'ClassyFire_level5', 'ClassyFire_most_specific')
  } else if (terms == 'custom'){
    term_levels = term_levels
  }
  for(term_level in term_levels){
    for(list in names(comp.list)){
      # Calculate enrichment score for all terms
      res <- CanopusLevelEnrichmentAnal(mmo = mmo, list_test = comp.list[[list]],sig=FALSE, pthr = pthr, representation = representation, term_level = term_level, pval = pval)
      res <- res |> mutate(comp = list)
      df.EA <- dplyr::bind_rows(df.EA, res)
      # get terms that are at least once enriched in one comparison
      res.sig <- CanopusLevelEnrichmentAnal(mmo = mmo, list_test = comp.list[[list]],sig=TRUE, pthr = pthr, representation = representation, term_level = term_level, pval = pval)
      sig.terms <- append(sig.terms, res.sig$term)
    }
    sig.terms <- unique(sig.terms)
    df.EA.sig <- df.EA |> filter(.data$term %in% sig.terms)
    if (pval == 'pval'){
      df.EA.sig <- df.EA.sig |>
        mutate(label = cut(
          .data$pval,
          breaks = c(0,0.001, 0.01, 0.05, 0.1, 1),
          labels = c("***", "**", "*", ".", "")
        ))
      plot <- ggplot(data = df.EA.sig, aes(x = .data$comp, y = .data$term, label = .data$label))+
        geom_point(aes(size = .data$subsetcount, color = .data$pval))+
        geom_text()+
        scale_size_area(name = 'Count', max_size = 10)+
        scale_color_gradient2(low = 'red', high = 'grey', mid = 'grey', midpoint = 0.4)+
        theme_minimal()+
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))+
        xlab('Comparisons')+
        ylab('Chemical classes')+
        facet_grid(term_level ~ ., scales = 'free_y', space = 'free', switch = 'y')
    } else if (pval == 'fdr'){
      df.EA.sig <- df.EA.sig |>
        mutate(label = cut(
          .data$fdr,
          breaks = c(0,0.001, 0.01, 0.05, 0.1, 1),
          labels = c("***", "**", "*", ".", "")
        ))
      plot <- ggplot(data = df.EA.sig, aes(x = .data$comp, y = .data$term, label = .data$label))+
        geom_point(aes(size = .data$subsetcount, color = .data$fdr))+
        geom_text()+
        scale_size_area(name = 'Count', max_size = 10)+
        scale_color_gradient2(low = 'red', high = 'grey', mid = 'grey', midpoint = 0.4)+
        theme_minimal()+
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))+
        xlab('Comparisons')+
        ylab('Chemical classes')+
        facet_grid(term_level ~ ., scales = 'free_y', space = 'free', switch = 'y')
    }
  }
  plot
  if (save_output){
    ggsave(paste0(outdir, '.pdf'), width = width, height = height)
    readr::write_csv(df.EA, paste0(outdir, '.csv'))
    readr::write_csv(df.EA.sig, paste0(outdir, '_sig.csv'))
  }
  return(list(plot = plot, df = df.EA))
}


#' Metabolite Set Enrichment Analysis (MSEA)
#'
#' This function performs Metabolite Set Enrichment Analysis (MSEA) using the fgsea package.
#' It takes a ranked list of feature scores and tests for enrichment of metabolite sets based on Canopus-predicted terms.
#' The results are saved as a CSV file and a PDF plot.
#' @param mmo The \code{mmo} with sirius annotation and normalized data
#' @param feature_id A vector of feature ids corresponding to the feature scores
#' @param feature_score A vector of feature scores (e.g., log2 fold changes)
#' @param term_level The level of term to use for enrichment analysis.
#'               Options are 'NPC_pathway', 'NPC_superclass', 'NPC_class',
#'              'ClassyFire_superclass', 'ClassyFire_class', 'ClassyFire_subclass',
#'             'ClassyFire_level5', or 'ClassyFire_most_specific' (default: 'NPC_class')
#' @param pthr The threshold for adjusted p-value to be considered significant (default: 0.05)
#' @param outdir The directory to save the output files (default: 'MSEA')
#' @param width The width of the output plot in inches (default: 8)
#' @param height The height of the output plot in inches (default: 12)
#' @param sig A logical value indicating whether to return only significant terms (default: FALSE)
#' @param save_output A logical value indicating whether to save the output plot (default: TRUE)
#' @return A list of the plot and the enrichment results
#' @export
#' @examplesIf FALSE
#' # Perform MSEA using NPC_class level
#' MSEA(
#'  mmo, feature_name = rownames(DE_results), feature_score = DE_results$log2FoldChange,
#'  term_level = 'NPC_class', pthr = 0.05, outdir = 'MSEA_NPC_class',
#'  width = 8, height = 12, sig = FALSE
#' )
MSEA <- function(mmo, feature_id, feature_score, term_level = 'NPC_class', pthr = 0.05, outdir = 'MSEA', width = 8, height = 12, sig = FALSE, save_output = TRUE){
  # Create a named vector of feature scores
  .require_pkg("fgsea")
  ranked_list <- feature_score
  names(ranked_list) <- feature_id
  ranked_list <- sort(ranked_list, decreasing = TRUE)

  # Retrieve metabolite sets based on the specified term level
  if(term_level == 'NPC_class'){
    metabolite_sets <- split(mmo$sirius_annot$id, mmo$sirius_annot[['NPC#class']])
  } else if (term_level == 'NPC_superclass'){
    metabolite_sets <- split(mmo$sirius_annot$id, mmo$sirius_annot[['NPC#superclass']])
  } else if (term_level == 'NPC_pathway'){
    metabolite_sets <- split(mmo$sirius_annot$id, mmo$sirius_annot[['NPC#pathway']])
  } else if (term_level == "ClassyFire_superclass") {
    metabolite_sets <- split(mmo$sirius_annot$id, mmo$sirius_annot[['ClassyFire#superclass']])
  } else if (term_level == "ClassyFire_class") {
    metabolite_sets <- split(mmo$sirius_annot$id, mmo$sirius_annot[['ClassyFire#class']])
  } else if (term_level == "ClassyFire_subclass") {
    metabolite_sets <- split(mmo$sirius_annot$id, mmo$sirius_annot[['ClassyFire#subclass']])
  } else if (term_level == "ClassyFire_level5") {
    metabolite_sets <- split(mmo$sirius_annot$id, mmo$sirius_annot[['ClassyFire#level 5']])
  } else if (term_level == "ClassyFire_most_specific") {
    metabolite_sets <- split(mmo$sirius_annot$id, mmo$sirius_annot[['ClassyFire#most specific class']])
  } else {
    stop("Invalid term level. Please choose a valid term level.")
  }
  msea_res <- fgsea::fgsea(pathways = metabolite_sets,
                       stats    = ranked_list,
                       minSize  = 5,   # minimum number of features in a class
                       maxSize  = 1500,
                       nPermSimple = 10000)
  msea_res <- msea_res |> arrange(.data$padj)
  # readr::write_csv(msea_res, paste0(outdir,'_', term_level,'_results.csv'))
  if (sig) {
    msea_res <- msea_res |> filter(.data$padj < pthr)
  }
  plot <- ggplot(msea_res, aes(x = reorder(.data$pathway, .data$NES), y = .data$NES)) +
    geom_point(shape = 21, aes(color = .data$padj < 0.05, size = .data$size, fill = -log(.data$padj)), stroke = 1) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    scale_fill_gradient(low = "grey", high = "red") +
    scale_color_manual(values = c("TRUE" = 'black', "FALSE" = 'white')) +
    guides(shape = "none") +
    labs(x = "Metabolite Class", y = "Normalized Enrichment Score (NES)", title = "MSEA Results", color = "-log10(padj)", size = "Set Size") +
    theme_classic() +
    theme(legend.position = "top", axis.text.y = element_text(size = 6))
  plot
  if (save_output){
    ggsave(paste0(outdir,'_', term_level,'.pdf'), width = width, height = height)
    readr::write_csv(msea_res, paste0(outdir,'_', term_level,'_results.csv'))
  }
  return (list(plot = plot, df = msea_res))
}


#' FeaturePhenotypeCorrelation
#'
#' This function performs correlation analysis of a specific feature against a phenotype in the metadata.
#' It can use linear mixed models (LMM), simple linear regression (LM), or Pearson correlation.
#' The default regression line of the plot uses linear model
#'
#' @param mmo The \code{mmo} with feature data and metadata
#' @param feature_id The id of the feature to analyze
#' @param phenotype The name of the phenotype  in the metadata
#' @param groups A vector of group names from the metadata containing phenotype data
#' @param model The type of regression model to use. Options are 'lmm' for linear mixed model, 'lm' for simple linear regression, or 'pearson' for Pearson correlation (default: 'lmm')
#' @param normalization The normalization method to use for feature data. Options are 'None', 'Log', 'Meancentered', or 'Z' (default: 'Z')
#' @param outdir The directory to save the output files
#' @param width The width of the output plot in inches (default: 6)
#' @param height The height of the output plot in inches (default: 6)
#' @param save_output A logical value indicating whether to save the output plot (default: TRUE)
#' @return A list of the plot and the raw data
#' @export

FeaturePhenotypeCorrelation <- function(mmo, feature_id, phenotype, groups, model = 'lm', normalization = 'Z', outdir = 'FeaturePhenotypeCorrelation', width = 6, height = 6, save_output = TRUE){
  .require_pkg("ggrepel")
  feature <- GetNormFeature(mmo, normalization)
  metadata <- mmo$metadata

  # Get phenotype phenotype from the metadata, get the feature value from the feature matrix, then combine
  phenotype.df <- data.frame(sample = metadata$sample, group = metadata$group, phenotype = metadata[,phenotype]) |> filter(.data$group %in% groups)
  feature_df <- data.frame(sample = colnames(feature[,-(1:2)]), feature_value = as.numeric(feature[feature$id == feature_id, -(1:2)]))
  combined_df <- merge(phenotype.df, feature_df, by='sample')

  # Perform linear mixed model or simple linear regression
  if (model == 'lmm'){
    fit <- lme4::lmer(combined_df$phenotype ~ combined_df$feature_value + (1|combined_df$group))
    p_value <- summary(fit)$coefficients[2, 5]
  } else if (model == 'lm'){
    fit <- lm(combined_df$phenotype ~ combined_df$feature_value)
    p_value <- summary(fit)$coefficients[2, 4]
  } else if (model %in% c('pearson', 'spearman', 'kendall')){
    correlation <- cor.test(combined_df$phenotype, combined_df$feature_value, method = model)
    p_value <- correlation$p.value
  } else {
    stop("Invalid model type. Please use 'lmm', 'lm', 'pearson', 'spearman' or 'kendall'")
  }

  # Plot the fit using ggplot
  plot <-ggplot(combined_df, aes(x = .data$feature_value, y = .data$phenotype, color = .data$group)) +
    geom_point(size = 3) +
    geom_smooth(method = "lm", se = TRUE, color = "black") +
    ggrepel::geom_text_repel(aes(label = sample), size = 2.5, show.legend = FALSE) +
    theme_classic() +
    labs(title = paste("Regression of", feature_id, "against", phenotype, "phenotype"),
         x = "Feature Value",
         y = "phenotype") +
    theme(legend.position = "right") +
    annotate("text", x = Inf, y = Inf, label = paste("p-value:", signif(p_value, digits = 4)),
             hjust = 1.1, vjust = 1.1, size = 3, color = "black")
  plot
  if (save_output){
    ggsave(paste0(outdir,'_', feature_id,'_vs_', phenotype,'.pdf'), height = height, width = width)
    readr::write_csv(combined_df, paste0(outdir,'_', feature_id,'_vs_', phenotype,'_data.csv'))
  }
  return (list(plot = plot, df = combined_df))
}


#' Screen feature-phenotype correlation
#'
#' Use metadata-provided variables (any phenotypes or environmental variables) to screen feature-phenotype correlation
#' linear model, linear mixed model (using groups as random effect), or correlation (Pearson, Spearman, Kendall) are supported
#'
#' @param mmo The \code{mmo} with feature data and metadata
#' @param phenotype The name of the phenotype in the metadata
#' @param groups A vector of group names from the metadata containing phenotype data
#' @param model The type of regression model to use. Options are 'lmm' for linear mixed model, 'lm' for simple linear regression, or 'pearson', 'spearman', 'kendall' for correlation (default: 'lm')
#' @param normalization The normalization method to use for feature data. Options are 'None', 'Log', 'Meancentered', or 'Z' (default: 'Z')
#' @return A list of the plot and the raw data
#' @export
ScreenFeaturePhenotypeCorrelation <- function(mmo, phenotype, groups, model = 'lm', normalization = 'None'){
  feature  <- GetNormFeature(mmo, normalization)
  metadata <- mmo$metadata

  if (missing(groups)) {
    groups <- unique(metadata$group)
    message("'groups' is not provided, using all groups")
  }
  phenotype_df <- data.frame(
    sample   = metadata$sample,
    group    = metadata$group,
    phenotype = metadata[, phenotype]
  ) |> filter(.data$group %in% groups)

  # Align feature matrix columns to phenotype samples
  samp_order <- intersect(colnames(feature)[-(1:2)], phenotype_df$sample)
  feat_mat   <- as.matrix(feature[, samp_order, drop = FALSE])
  y          <- phenotype_df$phenotype[match(samp_order, phenotype_df$sample)]

  # Samples with no phenotype value carry no information for any feature.
  # cor.test() and lm() drop them per feature; dropping them once up front is
  # equivalent and keeps them out of the matrix product below, which has no
  # na.rm and would otherwise turn every feature's covariance into NA.
  pheno_ok <- !is.na(y)
  if (!any(pheno_ok)) {
    stop("No samples have a non-missing '", phenotype, "' value in the selected groups.", call. = FALSE)
  }
  if (!all(pheno_ok)) {
    message(sum(!pheno_ok), " sample(s) dropped: no '", phenotype, "' value recorded.")
    feat_mat   <- feat_mat[, pheno_ok, drop = FALSE]
    samp_order <- samp_order[pheno_ok]
    y          <- y[pheno_ok]
  }
  n_feat     <- nrow(feat_mat)
  n          <- length(y)

  # Pre-allocate result vectors (avoids O(n^2) rbind accumulation)
  coeff_vec <- numeric(n_feat)
  pval_vec  <- numeric(n_feat)

  if (model %in% c('pearson', 'lm')) {
    # Vectorised correlation across all features at once. Features carrying a
    # missing value are handled one at a time below, because the pairwise
    # deletion they need changes both n and sd(y) for that feature alone.
    X      <- feat_mat
    y_use  <- y
    y_c    <- y_use - mean(y_use)
    X_c    <- sweep(X, 1L, rowMeans(X, na.rm = TRUE), `-`)
    n_obs  <- rowSums(!is.na(X))
    cov_xy <- as.numeric(X_c %*% y_c) / pmax(n_obs - 1L, 1L)
    sd_x   <- sqrt(rowSums(X_c^2L, na.rm = TRUE) / pmax(n_obs - 1L, 1L))
    sd_y   <- sd(y_use)
    r      <- cov_xy / (sd_x * sd_y)
    # A zero-variance feature gives 0/0. cor.test() reports NA for those, so
    # keep NA rather than reporting a correlation of exactly zero, which would
    # be indistinguishable from a real result.
    r[is.nan(r)] <- NA_real_

    coeff_vec <- if (model == 'lm') r * sd_y / pmax(sd_x, .Machine$double.eps) else r
    t_stat    <- r * sqrt(pmax(n - 2L, 1L) / pmax(1 - r^2, .Machine$double.eps))
    pval_vec  <- 2 * pt(-abs(t_stat), df = pmax(n - 2L, 1L))

    incomplete <- which(n_obs < n)
    for (i in incomplete) {
      xi <- feat_mat[i, ]
      ok <- !is.na(xi)
      if (sum(ok) < 3L || stats::sd(xi[ok]) == 0) {
        coeff_vec[i] <- NA_real_; pval_vec[i] <- NA_real_; next
      }
      res <- stats::cor.test(y[ok], xi[ok], method = "pearson")
      coeff_vec[i] <- if (model == 'lm') {
        unname(res$estimate) * stats::sd(y[ok]) / stats::sd(xi[ok])
      } else unname(res$estimate)
      pval_vec[i] <- res$p.value
    }

  } else if (model == 'spearman') {
    # cor.test() gives the exact (AS 89) p-value for small samples; the normal
    # approximation used by a vectorised rank correlation disagrees with it
    # enough to flip significance calls, so go through cor.test() per feature.
    for (i in seq_len(n_feat)) {
      xi <- feat_mat[i, ]
      ok <- !is.na(xi)
      if (sum(ok) < 3L || stats::sd(xi[ok]) == 0) {
        coeff_vec[i] <- NA_real_; pval_vec[i] <- NA_real_; next
      }
      res <- suppressWarnings(stats::cor.test(y[ok], xi[ok], method = "spearman"))
      coeff_vec[i] <- unname(res$estimate)
      pval_vec[i]  <- res$p.value
    }

  } else if (model %in% c('lmm', 'kendall')) {
    # Per-feature loop with pre-allocated vectors (no rbind)
    for (i in seq_len(n_feat)) {
      combined_df <- data.frame(
        sample        = samp_order,
        phenotype     = y,
        feature_value = feat_mat[i, ],
        group         = phenotype_df$group[match(samp_order, phenotype_df$sample)]
      )
      if (model == 'lmm') {
        fit           <- lme4::lmer(phenotype ~ feature_value + (1|group), data = combined_df)
        pval_vec[i]   <- summary(fit)$coefficients[2, 5]
        coeff_vec[i]  <- lme4::fixef(fit)[2]
      } else {
        res           <- cor.test(combined_df$phenotype, combined_df$feature_value, method = 'kendall')
        pval_vec[i]   <- res[[3]]
        coeff_vec[i]  <- res[[4]]
      }
    }
  } else {
    stop("Invalid model type. Use 'lmm', 'lm', 'pearson', 'spearman', or 'kendall'.")
  }

  data.frame(feature = feature$id, coefficient = coeff_vec, p_value = pval_vec,
             stringsAsFactors = FALSE)
}


#' GetPerformanceFeatureRegression
#'
#' This function performs linear regression analysis of all features against a phenotype performance in the metadata.
#' @param mmo The \code{mmo} with feature data and metadata
#' @param phenotype The name of the phenotype performance in the metadata
#' @param groups A vector of group names from the metadata containing performance data
#' @param DAM.list A list of DAMs to tag features
#' @param comparisons A list of pairwise comparisons to add fold change columns
#' @param normalization The normalization method to use for feature data. Options are 'None', 'Log', 'Meancentered', or 'Z' (default: 'None')
#' @export
#' @return A data frame containing regression results for each feature, including effect size, p-value, and fold change columns for specified comparisons.
#'
GetPerformanceFeatureRegression <- function(mmo, phenotype, groups, DAM.list, comparisons, normalization = 'None'){
  feature     <- GetNormFeature(mmo, normalization)
  metadata    <- mmo$metadata
  phenotype.df <- data.frame(
    sample = metadata$sample, group = metadata$group,
    performance = metadata[, phenotype]
  ) |> filter(.data$group %in% groups)

  # Align feature matrix to phenotype samples
  samp_order <- intersect(colnames(feature)[-(1:2)], phenotype.df$sample)
  feat_mat   <- as.matrix(feature[, samp_order, drop = FALSE])
  y          <- phenotype.df$performance[match(samp_order, phenotype.df$sample)]

  # Drop samples with no phenotype value before the matrix product below, which
  # has no na.rm and would otherwise propagate one missing value into every
  # feature's covariance. lm() discards these rows per feature anyway.
  pheno_ok <- !is.na(y)
  if (!any(pheno_ok)) {
    stop("No samples have a non-missing '", phenotype, "' value in the selected groups.", call. = FALSE)
  }
  if (!all(pheno_ok)) {
    message(sum(!pheno_ok), " sample(s) dropped: no '", phenotype, "' value recorded.")
    feat_mat   <- feat_mat[, pheno_ok, drop = FALSE]
    samp_order <- samp_order[pheno_ok]
    y          <- y[pheno_ok]
  }
  n_feat     <- nrow(feat_mat)
  n          <- length(y)

  # Vectorised simple linear regression via Pearson correlation
  y_c      <- y - mean(y)
  X_c      <- sweep(feat_mat, 1L, rowMeans(feat_mat, na.rm = TRUE), `-`)
  n_obs    <- rowSums(!is.na(feat_mat))
  cov_xy   <- as.numeric(X_c %*% y_c) / pmax(n_obs - 1L, 1L)
  var_x    <- rowSums(X_c^2L, na.rm = TRUE) / pmax(n_obs - 1L, 1L)
  eff_size <- cov_xy / pmax(var_x, .Machine$double.eps)   # lm slope = cov/var_x
  r        <- cov_xy / (sqrt(var_x) * sd(y))
  t_stat   <- r * sqrt(pmax(n - 2L, 1L) / pmax(1 - r^2, .Machine$double.eps))
  p_value  <- 2 * pt(-abs(t_stat), df = pmax(n - 2L, 1L))
  # A zero-variance feature has no slope to report; keep NA in both columns
  # rather than pairing an NA estimate with a confident-looking p of 1.
  degenerate <- is.nan(r) | is.na(r)
  p_value[is.nan(t_stat)] <- NA_real_
  eff_size[degenerate] <- NA_real_
  p_value[degenerate]  <- NA_real_

  # Features carrying a missing value need pairwise deletion, which changes
  # both n and sd(y) for that feature alone; handle those one at a time.
  for (i in which(n_obs < n)) {
    xi <- feat_mat[i, ]
    ok <- !is.na(xi)
    if (sum(ok) < 3L || stats::sd(xi[ok]) == 0) {
      eff_size[i] <- NA_real_; p_value[i] <- NA_real_; next
    }
    res <- stats::cor.test(y[ok], xi[ok], method = "pearson")
    eff_size[i] <- unname(res$estimate) * stats::sd(y[ok]) / stats::sd(xi[ok])
    p_value[i]  <- res$p.value
  }

  # DAM tags — vectorial lookup across all lists
  tags <- rep("else", n_feat)
  for (list_name in names(DAM.list)) {
    hits <- which(feature$id %in% DAM.list[[list_name]])
    tags[hits] <- list_name
  }

  regression_results <- data.frame(
    feature     = feature$id,
    effect.size = eff_size,
    p_value     = p_value,
    tag         = tags,
    stringsAsFactors = FALSE
  )
  for (comparison in comparisons) {
    fc_column <- paste(comparison, "log2FC", sep = "_")
    regression_results[[fc_column]] <- mmo$pairwise[[fc_column]][
      match(regression_results$feature, mmo$pairwise$feature)]
  }
  return(regression_results)
}

#' GetPerformanceFeatureLMM
#'
#' This function performs linear mixed model analysis of all features against a phenotype performance in the metadata.
#'
#' @param mmo The \code{mmo} with feature data and metadata
#' @param phenotype The name of the phenotype performance in the metadata
#' @param groups A vector of group names from the metadata containing performance data
#' @param DAM.list A list of DAMs to tag features
#' @param comparisons A list of pairwise comparisons to add fold change columns
#' @param normalization The normalization method to use for feature data. Options are 'None', 'Log', 'Meancentered', or 'Z' (default: 'Z')
#' @return A data frame containing regression results for each feature, including effect size, p-value, and fold change columns for specified comparisons.
#' @export
GetPerformanceFeatureLMM <- function(mmo, phenotype, groups, DAM.list, comparisons, normalization = 'Z'){
  feature      <- GetNormFeature(mmo, normalization)
  metadata     <- mmo$metadata
  phenotype.df <- data.frame(
    sample = metadata$sample, group = metadata$group,
    performance = metadata[, phenotype]
  ) |> filter(.data$group %in% groups)

  samp_order <- intersect(colnames(feature)[-(1:2)], phenotype.df$sample)
  feat_mat   <- as.matrix(feature[, samp_order, drop = FALSE])
  n_feat     <- nrow(feat_mat)

  # Pre-allocate result vectors (avoids O(n^2) rbind accumulation)
  eff_size_vec <- numeric(n_feat)
  p_value_vec  <- numeric(n_feat)

  grp_col <- phenotype.df$group[match(samp_order, phenotype.df$sample)]
  y_col   <- phenotype.df$performance[match(samp_order, phenotype.df$sample)]

  for (i in seq_len(n_feat)) {
    combined_df     <- data.frame(performance  = y_col,
                                  feature_value = feat_mat[i, ],
                                  group         = grp_col)
    lmm_fit          <- lme4::lmer(performance ~ feature_value + (1|group), data = combined_df)
    eff_size_vec[i]  <- lme4::fixef(lmm_fit)[2]
    p_value_vec[i]   <- summary(lmm_fit)$coefficients[2, 5]
  }

  # DAM tags — vectorial lookup
  tags <- rep("else", n_feat)
  for (list_name in names(DAM.list)) {
    hits        <- which(feature$id %in% DAM.list[[list_name]])
    tags[hits]  <- list_name
  }

  regression_results <- data.frame(
    feature     = feature$id,
    effect.size = eff_size_vec,
    p_value     = p_value_vec,
    tag         = tags,
    stringsAsFactors = FALSE
  )
  for (comparison in comparisons) {
    fc_column <- paste(comparison, "log2FC", sep = "_")
    regression_results[[fc_column]] <- mmo$pairwise[[fc_column]][
      match(regression_results$feature, mmo$pairwise$feature)]
  }
  return(regression_results)
}

#' GetPerformanceFeatureCorrelation
#'
#'
#' This function calculates the Pearson correlation between each feature and a specified phenotype performance in the metadata.
#' @param mmo The \code{mmo} with feature data and metadata
#' @param phenotype The name of the phenotype performance in the metadata
#' @param groups A vector of group names from the metadata containing performance data
#' @param DAM.list A list of DAMs to tag features
#' @param comparisons A list of pairwise comparisons to add fold change columns
#' @param cor_method The correlation method to use. Options are 'pearson', 'spearman', or 'kendall' (default: 'pearson')
#' @param normalization The normalization method to use for feature data. Options are 'None', 'Log', 'Meancentered', or 'Z' (default: 'None')
#' @return A data frame containing correlation results for each feature, including effect size (correlation coefficient), p-value, and fold change columns for specified comparisons.
#' @export
GetPerformanceFeatureCorrelation <- function(mmo, phenotype, groups, DAM.list, comparisons, cor_method = 'pearson', normalization = 'None'){
  feature      <- GetNormFeature(mmo, normalization)
  metadata     <- mmo$metadata
  phenotype.df <- data.frame(
    sample = metadata$sample, group = metadata$group,
    performance = metadata[, phenotype]
  ) |> filter(.data$group %in% groups)

  samp_order <- intersect(colnames(feature)[-(1:2)], phenotype.df$sample)
  feat_mat   <- as.matrix(feature[, samp_order, drop = FALSE])
  y          <- phenotype.df$performance[match(samp_order, phenotype.df$sample)]

  # Drop samples with no phenotype value before the matrix product below, which
  # has no na.rm and would otherwise propagate one missing value into every
  # feature's covariance. cor.test() discards these pairs per feature anyway.
  pheno_ok <- !is.na(y)
  if (!any(pheno_ok)) {
    stop("No samples have a non-missing '", phenotype, "' value in the selected groups.", call. = FALSE)
  }
  if (!all(pheno_ok)) {
    message(sum(!pheno_ok), " sample(s) dropped: no '", phenotype, "' value recorded.")
    feat_mat   <- feat_mat[, pheno_ok, drop = FALSE]
    samp_order <- samp_order[pheno_ok]
    y          <- y[pheno_ok]
  }
  n_feat     <- nrow(feat_mat)
  n          <- length(y)

  # Pre-allocate result vectors
  eff_size_vec <- numeric(n_feat)
  p_value_vec  <- numeric(n_feat)

  if (cor_method == 'pearson') {
    X      <- feat_mat
    y_c    <- y - mean(y)
    X_c    <- sweep(X, 1L, rowMeans(X, na.rm = TRUE), `-`)
    n_obs  <- rowSums(!is.na(X))
    cov_xy <- as.numeric(X_c %*% y_c) / pmax(n_obs - 1L, 1L)
    sd_x   <- sqrt(rowSums(X_c^2L, na.rm = TRUE) / pmax(n_obs - 1L, 1L))
    sd_y   <- sd(y)
    r      <- cov_xy / (sd_x * sd_y)
    # A zero-variance feature gives 0/0; cor.test() reports NA for those.
    r[is.nan(r)] <- NA_real_
    t_stat       <- r * sqrt(pmax(n - 2L, 1L) / pmax(1 - r^2, .Machine$double.eps))
    eff_size_vec <- r
    p_value_vec  <- 2 * pt(-abs(t_stat), df = pmax(n - 2L, 1L))

    # Features carrying a missing value need pairwise deletion, which changes
    # both n and sd(y) for that feature alone; handle those one at a time.
    for (i in which(n_obs < n)) {
      xi <- feat_mat[i, ]; ok <- !is.na(xi)
      if (sum(ok) < 3L || stats::sd(xi[ok]) == 0) {
        eff_size_vec[i] <- NA_real_; p_value_vec[i] <- NA_real_; next
      }
      res <- stats::cor.test(y[ok], xi[ok], method = "pearson")
      eff_size_vec[i] <- unname(res$estimate); p_value_vec[i] <- res$p.value
    }
  } else {
    # spearman and kendall: cor.test() supplies the exact small-sample p-value,
    # which a rank-correlation t approximation does not reproduce.
    for (i in seq_len(n_feat)) {
      xi <- feat_mat[i, ]; ok <- !is.na(xi)
      if (sum(ok) < 3L || stats::sd(xi[ok]) == 0) {
        eff_size_vec[i] <- NA_real_; p_value_vec[i] <- NA_real_; next
      }
      res             <- suppressWarnings(stats::cor.test(y[ok], xi[ok], method = cor_method))
      eff_size_vec[i] <- unname(res$estimate)
      p_value_vec[i]  <- res$p.value
    }
  }

  # DAM tags — vectorial lookup
  tags <- rep("else", n_feat)
  for (list_name in names(DAM.list)) {
    hits       <- which(feature$id %in% DAM.list[[list_name]])
    tags[hits] <- list_name
  }

  regression_results <- data.frame(
    feature     = feature$id,
    effect.size = eff_size_vec,
    p_value     = p_value_vec,
    tag         = tags,
    stringsAsFactors = FALSE
  )
  for (comparison in comparisons) {
    fc_column <- paste(comparison, "log2FC", sep = "_")
    regression_results[[fc_column]] <- mmo$pairwise[[fc_column]][
      match(regression_results$feature, mmo$pairwise$feature)]
  }
  return(regression_results)
}


#' PlotFoldchangeResistanceRegression
#'
#' This function plots the regression results of a feature against a fold change in resistance, including
#' regression line, p-value, and R-squared value.
#'
#' @param performance_regression The regression results data frame containing effect size, fold change, and tag. The output from GetPerformanceFeatureRegression, GetPerformanceFeatureLMM, or GetPerformanceFeatureCorrelation.
#' @param fold_change The name of the fold change column in the performance_regression dataframe
#' @param color A vector of colors for the points in the plot
#' @param outdir The output file path for the regression plot
#' @param width The width of the output plot in inches (default: 6)
#' @param height The height of the output plot in inches (default: 6)
#' @param save_output A logical value indicating whether to save the output plot (default: TRUE)
#' @return A list containing the regression plot and the performance regression data frame
#' @export
PlotFoldchangeResistanceRegression <- function(performance_regression, fold_change, color, outdir, width = 6, height = 6, save_output = TRUE){
  ind_fit <- lm(data = performance_regression, formula = as.formula(paste("-effect.size ~", fold_change)))
  summary_fit <- summary(ind_fit)
  p_value <- summary_fit$coefficients[2, 4]
  r_squared <- summary_fit$r.squared

  plot <- ggplot(performance_regression, aes(x = !!rlang::sym(fold_change), y = -.data$effect.size)) +
    geom_point(size = 0.5, aes(color = .data$tag)) +
    geom_smooth(method = "lm", se = TRUE, color = "black", level = 0.95) +
    xlab(fold_change) +
    ylab('-effect.size') +
    scale_color_manual(values = color) +
    theme_classic() +
    annotate("text", x = Inf, y = Inf, label = paste("p-value:", round(p_value, 500), "\nR-squared:", round(r_squared, 4)),
            hjust = 1.1, vjust = 1.1, size = 3, color = "black")+
    geom_hline(yintercept = 0, linetype = "dashed", color = "black")
  plot
  if (save_output){
    ggsave(paste0(outdir, ".png"), plot, height = height, width = width)
    readr::write_csv(performance_regression, file = paste0(outdir, ".csv"))
  }
  return(list(plot = plot, df = performance_regression))
}

#' PlotFoldchangeResistanceRegression_t
#'
#' This function plots the regression results of a feature against a fold change in resistance, including
#' regression line, p-value, and R-squared value. Transposed version of PlotFoldchangeResistanceRegression.
#'
#' @param performance_regression The regression results data frame containing effect size, fold change, and tag. The output from GetPerformanceFeatureRegression, GetPerformanceFeatureLMM, or GetPerformanceFeatureCorrelation.
#' @param fold_change The name of the fold change column in the performance_regression dataframe
#' @param color A vector of colors for the points in the plot
#' @param output_dir The output file path for the regression plot
#' @param save_output A logical value indicating whether to save the output plot (default: TRUE)
#' @param width The width of the output plot in inches (default: 6)
#' @param height The height of the output plot in inches (default: 6)
#' @return A list containing the regression plot and the performance regression data frame
#' @export
PlotFoldchangeResistanceRegression_t <- function(performance_regression, fold_change, color, output_dir, save_output = TRUE, width = 6, height = 6){
  ind_fit <- lm(data = performance_regression, formula = as.formula(paste(fold_change, "~ -effect.size")))
  summary_fit <- summary(ind_fit)
  p_value <- summary_fit$coefficients[4]
  r_squared <- summary_fit$r.squared

  plot <- ggplot(performance_regression, aes(x = -.data$effect.size, y = !!rlang::sym(fold_change))) +
    geom_point(size = 0.5, aes(color = .data$tag)) +
    geom_smooth(method = "lm", se = TRUE, color = "black", level = 0.95) +
    xlab('-effect.size') +
    ylab(fold_change) +
    scale_color_manual(values = color) +
    theme_classic() +
    annotate("text", x = Inf, y = Inf, label = paste("p-value:", round(p_value, 500), "\nR-squared:", round(r_squared, 4)),
            hjust = 1.1, vjust = 1.1, size = 3, color = "black")+
    geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
    geom_vline(xintercept = 0, linetype = "dashed", color = "black")
  if (save_output){
    ggsave(paste0(output_dir, ".png"), plot, height = height, width = width)
    readr::write_csv(performance_regression, file = paste0(output_dir, ".csv"))
  }
  return(list(plot = plot, df = performance_regression))
}

#' PlotFoldchangeResistanceQuad
#'
#' This function plots the fold change resistance in a quadrant plot, categorizing points into quadrants based on their effect size and fold change.
#' It also performs a binomial test to assess the distribution of points across quadrants.
#' @param performance_regression The regression results data frame containing effect size, fold change, and tag. The output from GetPerformanceFeatureRegression, GetPerformanceFeatureLMM, or GetPerformanceFeatureCorrelation.
#' @param fold_change The name of the fold change column in the performance_regression dataframe
#' @param color A vector of colors for the points in the plot
#' @param output_dir The output file path for the quadrant plot
#' @param save_output A logical value indicating whether to save the output plot (default: TRUE)
#' @param width The width of the output plot in inches (default: 6)
#' @param height The height of the output plot in inches (default: 6)
#' @return A list containing the quadrant plot and the performance regression data frame
#' @export
PlotFoldchangeResistanceQuad <- function(performance_regression, fold_change, color, output_dir, save_output = TRUE, width = 6, height = 6){
  performance_regression <- performance_regression |>
  mutate(
    quadrant = dplyr::case_when(
      -effect.size > 0 & !!rlang::sym(fold_change) > 0 ~ "Q1",
      -effect.size < 0 & !!rlang::sym(fold_change) < 0 ~ "Q3",
      -effect.size < 0 & !!rlang::sym(fold_change) > 0 ~ "Q2",
      -effect.size > 0 & !!rlang::sym(fold_change) < 0 ~ "Q4",
      TRUE ~ "Edge"             # For points on axes
    )
  )
  q_counts <- table(performance_regression$quadrant)
  q13 <- sum(q_counts[c("Q1", "Q3")], na.rm = TRUE)
  q24 <- sum(q_counts[c("Q2", "Q4")], na.rm = TRUE)
  binom_test <- binom.test(q13, q13+q24, p = 0.5, alternative = "two.sided")

  plot <- ggplot(performance_regression, aes(x = -.data$effect.size, y = !!sym(fold_change))) +
    geom_point(size = 0.5, aes(color = .data$tag)) +
    xlab('-effect.size') +
    ylab(fold_change) +
    scale_color_manual(values = color) +
    theme_classic() +
    annotate("text", x = Inf, y = Inf, label = paste("p-value:", round(binom_test[[3]], 500)),
            hjust = 1.1, vjust = 1.1, size = 3, color = "black")+
    geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
    geom_vline(xintercept = 0, linetype = "dashed", color = "black")
  plot
  if (save_output){
    ggsave(paste0(output_dir, ".png"), plot, height = height, width = width)
    readr::write_csv(performance_regression, file = paste0(output_dir, ".csv"))
  }
  return(list(plot = plot, df = performance_regression))
}
################### Singlevariate analyses ###################

#' Generate barplots for each feature and perform ANOVA
#'
#' This function generates bar plots for a specified feature across different groups in the metadata, performing ANOVA and Tukey's HSD test for post-hoc analysis.
#'
#' @param mmo The \code{mmo} containing metadata and feature data
#' @param ID_list A list of feature IDs to analyze. Use FeatureToID() to convert feature names to IDs.
#' @param outdir The output directory to save the bar plots and ANOVA results
#' @param normalization The normalization method to use for feature data. Options are 'None', 'Log', 'Meancentered', or 'Z' (default: 'None')
#' @param filter_group A boolean indicating whether to filter the feature values by a specific group list (default: FALSE)
#' @param group_list A list of groups to filter the feature values by, if filter_group is TRUE (default: NULL)
#' @param save_output A logical value indicating whether to save the output plot (default: TRUE)
#' @param width The width of the output plot in inches (default: 6)
#' @param height The height of the output plot in inches (default: 6)
#' @param colors A vector of colors for the groups (default: NULL)
#' @param posthoc A character string indicating the post-hoc test to perform (default: NULL)
#' @export
#' @return A list containing the bar plot and the ANOVA results
#' @examplesIf FALSE
#' AnovaBarPlot(mmo, ID_list = c("ID_1", "ID_2"), outdir = "output_directory", normalization = 'Z')
#' AnovaBarPlot(
#'  mmo, ID_list = c("ID_1", "ID_2"), outdir = "output_directory", normalization = 'Z',
#'  filter_group = TRUE, group_list = c("Group1", "Group2")
#' )
AnovaBarPlot <- function(mmo, ID_list, outdir, normalization = 'None', filter_group = FALSE, group_list = NULL,
                          save_output = TRUE, width = 6, height = 6, colors = NULL, posthoc = NULL) {
  .require_pkg("ggbeeswarm")
  .require_pkg("cld")
  if (filter_group){
    mmo <- filter_mmo(mmo, group_list = group_list)
  }
  # Extract metadata and feature data
  metadata <- mmo$metadata
  feature_data <- GetNormFeature(mmo, normalization)

  # Iterate through each feature ID
  for (target_id in ID_list) {
    # Extract feature values and merge with metadata
    feature_values <- feature_data |>
      filter(.data$id == target_id) |>
      dplyr::select(-.data$id, -.data$feature) |>
      t() |>
      as.data.frame()
    colnames(feature_values) <- "value"
    feature_values$sample <- rownames(feature_values)
    feature_values <- merge(feature_values, metadata, by = "sample")
    # Perform ANOVA
    anova <- anova_tukey_dunnett(feature_values, 'value ~ group')

    Letters <- NULL
    if(!is.null(posthoc)){
        if(posthoc == "Tukey"){
          res.Letters <- multcompView::multcompLetters4(anova$aov_res, anova$tukey_res)
          Letters <- data.frame(x = names(res.Letters$group$Letters),
                              label = as.character(res.Letters$group$Letters))
        } else if(posthoc == "Dunnett"){
          temp <- cld::make_cld(anova$dunnett_res)
          Letters <- data.frame(x = temp$group, label = temp$cld)
        } else{
          message("Invalid post-hoc method (labelling skipped)")
          Letters <- NULL
        }
    }


    # Generate bar plot
    plot <- ggplot(feature_values, aes(x = .data$group, y = .data$value, fill = .data$group)) +
      geom_bar(stat = "summary", fun = "mean", position = "dodge") +
      geom_errorbar(stat = "summary", fun.data = "mean_se", position = position_dodge(width = 0.9), width = 0.2) +
      ggbeeswarm::geom_beeswarm() +
      theme_classic() +
      labs(title = paste("Feature: ", target_id, "; Post-hoc: ", posthoc, sep=""), x = "Group", y = "Value") +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

    if(!is.null(colors)){ ## BE: ! if manual color is exist
        plot <- plot +
                scale_fill_manual(values = colors)
    }
    if(!is.null(Letters)){ ## BE: ! if post-hoc is..
      plot <- plot +
              geom_text(data = Letters, aes(x = .data$x, y = max(feature_values$value)*1.1, label = .data$label),
              inherit.aes = FALSE, size = 3, color = "black")
    }
    plot
    if (save_output){
      ggsave(file.path(outdir, paste0(target_id, "_barplot.pdf")), plot = plot, width = width, height = height)
      write_anova(anova, outdir = paste0(outdir,'/', target_id, '_anova.csv'), way = 'oneway')
    }
  }
  return(list(plot = plot, anova = anova))
}

#' ExportFeaturesToCSV
#'
#' This function exports selected features, their annotations, and pairwise comparisons to a CSV file.
#'
#' @param mmo The \code{mmo} containing feature data, annotations, and pairwise comparisons
#' @param id_list A list of feature names to filter and export
#' @param normalization The normalization method to use for feature data.
#'        Options are 'None', 'Log', 'Meancentered', or 'Z' (default: 'None')
#' @param output_dir The output directory to save the CSV file
#' @export
#' @examplesIf FALSE
#' ExportFeaturesToCSV(mmo,
#'                     id_list = Glucosinolates,
#'                      normalization = 'Z',
#'                      output_dir = 'output.csv')
#' ExportFeaturesToCSV(mmo, id_list = DAMs_up$control_vs_treatment1.up,
#'                          normalization = 'None',
#'                           output_dir = 'output.csv')
#'
ExportFeaturesToCSV <- function(mmo, id_list, normalization = 'None', output_dir){
  feature <- GetNormFeature(mmo, normalization = normalization) # Get normalized feature data
  # Filter the feature data, annotation, and DA analysis for the list provided
  selected_feature <- feature |> filter(.data$id %in% id_list)
  selected_pairwise <- mmo$pairwise |> filter(.data$id %in% id_list)
  # Merge all
  merged_df <- merge(mmo$sirius_annot, selected_feature, by = 'id')
  merged_df <- merge(merged_df, selected_pairwise, by = 'id')

  readr::write_csv(merged_df, output_dir)
}

#' GetRichness
#'
#' Sample-level richness: number of features present in each sample.
#' A feature is present if value > threshold; 0 never counts as present.
#'
#' @param feature_data Feature table with columns: id, feature, then sample columns
#' @param metadata Metadata table with sample and group columns
#' @param threshold Numeric; detection threshold for presence (default: 0)
#' @export
#' @return data.frame with columns: sample, group, richness
GetRichness <- function(
    feature_data,
    metadata,
    threshold = 0
) {
  sample_cols <- colnames(feature_data)[-c(1, 2)]
  x <- feature_data[, sample_cols, drop = FALSE]

  richness <- colSums((x > threshold) & !is.na(x))

  group_map <- metadata$group[match(sample_cols, metadata$sample)]

  data.frame(
    sample = sample_cols,
    group = group_map,
    richness = as.numeric(richness),
    stringsAsFactors = FALSE
  )
}


# DEPRECATED: superseded by GetFunctionalHillNumber().
# GetFunctionalHillNumber_derep() takes a dissimilarity matrix (D) and computes
# Rao Q as p'Dp. GetFunctionalHillNumber() takes a similarity matrix (S) and
# uses the algebraically equivalent p'Dp = 1 - p'Sp (valid when sum(p) = 1).
#' GetFunctionalHillNumber
#'
#' This function calculates the functional Hill number for a given \code{mmo}, normalization method, and distance metric.
#' See https://nph.onlinelibrary.wiley.com/doi/full/10.1111/nph.18685 for details of the functional Hill number calculation.
#'
#' @param feature Feature table with columns: id, feature, then sample columns
#' @param metadata Metadata table with sample and group columns
#' @param distance_matrix Feature distance matrix
#' @param q The order of the Hill number to calculate (default: 1).
#'        Larger q values give more weight to evenness portion of the hill number over richness.
#'        q = 0 treats abundance matrix as a presence absence matrix
#' @param threshold Numeric; detection threshold for presence (default: 0)
#' @param scale_dissim Boolean; whether to scale the distance matrix to be between 0 and 1 (default: TRUE)
#' @export
#' @return A data frame containing the functional Hill number for each group in the metadata, with columns for group and hill number.
GetFunctionalHillNumber_derep <- function(
    feature,
    metadata,
    distance_matrix,
    q = 1,
    threshold = 0,
    scale_dissim = TRUE
){
  # Treat values <= threshold (or NA) as absent
  feature_thr <- feature
  x_thr <- as.matrix(feature_thr[, -(1:2), drop = FALSE])
  x_thr[is.na(x_thr) | x_thr <= threshold] <- 0
  feature_thr[, -(1:2)] <- x_thr

  # Scale the  distance matrix to be between 0 and 1
  if(scale_dissim){
    scaled_dissimilarity <- distance_matrix / max(distance_matrix)
  } else {
    scaled_dissimilarity <- distance_matrix
  }
  # Calculate the relative proportions of each feature and reorder them to match the order of the distance matrix
  q.feature <- feature_thr |> filter(.data$id %in% colnames(scaled_dissimilarity))
  relative_proportions <- apply(q.feature[, -(1:2)], 2, function(x) {
    s <- sum(x, na.rm = TRUE)
    if (s <= 0) {
      rep(0, length(x))
    } else {
      x / s
    }
  })
  # apply() drops dimensions when there is only 1 column (e.g. pooled rarefaction);
  # restore as a named matrix so colnames() is never NULL downstream.
  if (is.null(dim(relative_proportions))) {
    nm <- colnames(q.feature[, -(1:2), drop = FALSE])
    relative_proportions <- matrix(relative_proportions, ncol = 1,
                                   dimnames = list(names(relative_proportions), nm))
  }
  rownames(relative_proportions) <- q.feature$id
  relative_proportions <- relative_proportions[rownames(scaled_dissimilarity), , drop = FALSE]
  scaled_dissimilarity <- as.matrix(scaled_dissimilarity)
  raoQ <- colSums(relative_proportions * (scaled_dissimilarity %*% relative_proportions))
  # Calculate Hill
  functional_hill_number <- c()
  if (q == 1){
    mask <- relative_proportions > 0
    Plog <- ifelse(mask, relative_proportions/raoQ * log(relative_proportions/raoQ), 0)
    DP <- scaled_dissimilarity %*% relative_proportions
    vals <- 2 * colSums(Plog * DP)
    functional_hill_number <- exp(-vals)
  } else {
    Pq <- (relative_proportions/raoQ)^q
    DPq <- scaled_dissimilarity %*% Pq
    vals <- colSums(Pq*DPq)
    functional_hill_number <- vals^(1/(1-q))
  }
  sample_names <- colnames(relative_proportions)
  names(functional_hill_number) <- sample_names
  # Get the group information
  groups <- metadata$group[match(sample_names, metadata$sample)]
  hill_df <- data.frame(
    sample = sample_names,
    group = groups,
    hill_number = as.numeric(functional_hill_number),
    value = as.numeric(functional_hill_number),
    stringsAsFactors = FALSE
  )
  return(hill_df)
}


#' GetFunctionalHillNumber (similarity-based)
#'
#' Calculates the functional Hill number using a pairwise feature \strong{similarity}
#' matrix (values 0--1, higher = more similar) from \code{GetSimMat()}.
#' Algebraically equivalent to \code{GetFunctionalHillNumber_derep()} but avoids
#' ever constructing the full dense distance matrix, making it safe for sparse
#' large-dataset similarity matrices. Key identity used: \eqn{p'Dp = 1 - p'Sp}
#' (exact when \eqn{\sum p_i = 1} and \eqn{D = 1 - S}).
#'
#' @param feature Feature table with columns: id, feature, then sample columns
#' @param metadata Metadata table with sample and group columns
#' @param sim_matrix Feature similarity matrix (dense or sparse dgCMatrix)
#' @param q The order of the Hill number (default: 1)
#' @param threshold Numeric; detection threshold for presence (default: 0)
#' @param scale_dissim Kept for interface compatibility; max(sim)=1 so scaling
#'   is a no-op, but passing FALSE still skips it
#' @export
#' @return A data frame with sample, group, hill_number, and value columns.
GetFunctionalHillNumber <- function(
    feature,
    metadata,
    sim_matrix,
    q = 1,
    threshold = 0,
    scale_dissim = TRUE
) {
  # Treat values <= threshold (or NA) as absent
  feature_thr <- feature
  x_thr <- as.matrix(feature_thr[, -(1:2), drop = FALSE])
  x_thr[is.na(x_thr) | x_thr <= threshold] <- 0
  feature_thr[, -(1:2)] <- x_thr

  # Scale similarity to [0,1]; max(sim)=1 (diagonal=self-similarity), so trivial
  if (scale_dissim) {
    scaled_similarity <- sim_matrix / max(sim_matrix)
  } else {
    scaled_similarity <- sim_matrix
  }

  # Calculate relative proportions and align to similarity matrix row order
  q.feature <- feature_thr |> filter(.data$id %in% colnames(scaled_similarity))
  relative_proportions <- apply(q.feature[, -(1:2)], 2, function(x) {
    s <- sum(x, na.rm = TRUE)
    if (s <= 0) rep(0, length(x)) else x / s
  })
  # apply() drops dimensions when there is only 1 column (e.g. pooled rarefaction);
  # restore as a named matrix so colnames() is never NULL downstream.
  if (is.null(dim(relative_proportions))) {
    nm <- colnames(q.feature[, -(1:2), drop = FALSE])
    relative_proportions <- matrix(relative_proportions, ncol = 1,
                                   dimnames = list(names(relative_proportions), nm))
  }
  rownames(relative_proportions) <- q.feature$id
  relative_proportions <- relative_proportions[rownames(scaled_similarity), , drop = FALSE]

  # Rao Q via similarity: p'Dp = 1 - p'Sp (exact identity when sum(p)=1, D=1-S).
  # No as.matrix() call here — sparse S passes directly into %*% without densification.
  SP   <- scaled_similarity %*% relative_proportions
  raoQ <- 1 - colSums(relative_proportions * SP)

  # Calculate Hill numbers
  functional_hill_number <- c()
  if (q == 1) {
    mask <- relative_proportions > 0
    Plog <- ifelse(mask, relative_proportions / raoQ * log(relative_proportions / raoQ), 0)
    # Dp = 1 - Sp  (valid because sum(p) = 1, so J%*%p = 1)
    DP   <- 1 - SP
    vals <- 2 * colSums(Plog * DP)
    functional_hill_number <- exp(-vals)
  } else {
    Pq  <- (relative_proportions / raoQ)^q
    SPq <- scaled_similarity %*% Pq
    # Dpq = colSums(Pq)*1 - Spq  (general: Pq does not sum to 1, so J%*%Pq != 1)
    DPq <- sweep(SPq, 2, colSums(Pq), function(sp, cs) cs - sp)
    vals <- colSums(Pq * DPq)
    functional_hill_number <- vals^(1 / (1 - q))
  }
  sample_names <- colnames(relative_proportions)
  names(functional_hill_number) <- sample_names
  # Get the group information
  groups <- metadata$group[match(sample_names, metadata$sample)]
  hill_df <- data.frame(
    sample = sample_names,
    group = groups,
    hill_number = as.numeric(functional_hill_number),
    value = as.numeric(functional_hill_number),
    stringsAsFactors = FALSE
  )
  return(hill_df)
}

#' GetHillNumbers
#'
#' This function calculates the Hill numbers for a given \code{mmo}, normalization method, and order of the Hill number without considering feature dissimilarity.
#'
#'
#' @param feature Feature table with columns: id, feature, then sample columns
#' @param metadata Metadata table with sample and group columns
#' @param q The order of the Hill number to calculate (default: 0)
#' @param threshold Numeric; detection threshold for presence (default: 0)
#' @export
#' @return A data frame containing the Hill number for each group in the metadata, with columns for group and hill number.
GetHillNumbers <- function(
    feature,
    metadata,
    q = 0,
    threshold = 0
) {
  x_thr <- as.matrix(feature[, -(1:2), drop = FALSE])
  x_thr[is.na(x_thr) | x_thr <= threshold] <- 0

  # All-zero samples produce NA diversity
  sample_sums <- colSums(x_thr, na.rm = TRUE)

  hill_numbers <- apply(feature[, -(1:2)], 2, function(x) {
    x <- ifelse(is.na(x) | x <= threshold, 0, x)
    if (sum(x, na.rm = TRUE) <= 0) {
      return(NA_real_)
    }
    p <- x / sum(x, na.rm = TRUE)
    if (q == 0) {
      return(sum(x > 0, na.rm = TRUE))
    } else if (q == 1) {
      p <- p[p > 0]
      return(exp(-sum(p * log(p), na.rm = TRUE)))
    } else {
      return((sum(p^q, na.rm = TRUE))^(1 / (1 - q)))
    }
  })
  hill_numbers[sample_sums <= 0] <- NA_real_

  sample_names <- names(hill_numbers)
  groups <- metadata$group[match(sample_names, metadata$sample)]
  hill_df <- data.frame(
    sample = sample_names,
    group = groups,
    hill_number = as.numeric(hill_numbers),
    value = as.numeric(hill_numbers),
    stringsAsFactors = FALSE
  )


  return(hill_df)
}

#' GetFaithPD
#'
#' This function calculates the Faith's phylogenetic diversity for a given \code{mmo} and distance metric,
#' to calculate chemically-informed richness
#'
#' @param feature Feature table with columns: id, feature, then sample columns
#' @param metadata Metadata table with sample and group columns
#' @param distance_matrix Feature distance matrix
#' @param threshold Numeric; detection threshold for presence (default: 0)
#' @export
#' @return A data frame containing the Faith's phylogenetic diversity for each group in the metadata, with columns for group and PD.
# DEPRECATED: superseded by GetFaithPD().
# GetFaithPD_derep() takes a dissimilarity matrix; GetFaithPD() takes a
# similarity matrix and adds a size guard before the unavoidable densification.
GetFaithPD_derep <- function(feature, metadata, distance_matrix, threshold = 0){
  .require_pkg("picante")
  .require_pkg("ape")
  ids <- feature$id
  x <- as.matrix(feature[, -(1:2), drop = FALSE])
  # Faith PD uses presence/absence; threshold defines "present"
  x <- ifelse(!is.na(x) & x > threshold, 1, 0)
  feature_t <- t(x)
  colnames(feature_t) <- ids
  tree <- ape::as.phylo(hclust(as.dist(distance_matrix), method = 'average'))
  pd_result <- picante::pd(feature_t, tree)
  pd_result$sample <- rownames(pd_result)
  pd_result$group <- metadata$group[match(pd_result$sample, metadata$sample)]
  pd_result$value <- pd_result$PD
  rownames(pd_result) <- NULL
  return(pd_result)
}


# sparse_single_phylo: builds a phylo tree from a sparse OR dense similarity
# matrix using the MST = single-linkage equivalence. Avoids O(n^2)
# densification — only explicit (non-zero) similarity edges are loaded into
# igraph; absent entries are treated as maximally distant (distance = M).
# Used by GetFaithPD and GetBetaDiversity when use_mst = TRUE or n > 10000.
sparse_single_phylo <- function(S, M = 1) {
  .require_pkg("igraph")
  .require_pkg("ape")

  S  <- methods::as(S, "TsparseMatrix")        # i, j, x triplets (0-based)
  up <- S@i < S@j                               # upper triangle only; excludes diagonal
  ii <- S@i[up] + 1L
  jj <- S@j[up] + 1L
  w  <- M - S@x[up]                            # convert similarity to distance
  n  <- nrow(S)

  g <- igraph::make_empty_graph(n = n, directed = FALSE)
  g <- igraph::add_edges(g, as.vector(rbind(ii, jj)))
  igraph::E(g)$weight <- w
  mst  <- igraph::mst(g)                        # C backend; spanning forest if disconnected
  # as_edgelist returns double indices; coerce once so parent[] stays integer throughout
  el   <- matrix(as.integer(igraph::as_edgelist(mst, names = FALSE)), ncol = 2L)
  ew   <- igraph::E(mst)$weight
  o    <- order(ew)
  el   <- el[o, , drop = FALSE]
  ew   <- ew[o]

  # union-find -> hclust(merge, height)
  parent <- seq_len(n)
  csize  <- rep(1L, n)
  cid    <- -seq_len(n)          # singletons are negative in hclust convention

  find <- function(x) {
    while (parent[x] != x) { parent[x] <<- parent[parent[x]]; x <- parent[x] }
    x
  }

  merge  <- matrix(0L, n - 1L, 2L)
  height <- numeric(n - 1L)
  k      <- 0L

  link <- function(a, b, h) {
    ra <- find(a); rb <- find(b)
    k <<- k + 1L
    merge[k, ] <<- c(cid[ra], cid[rb])
    height[k]  <<- h
    if (csize[ra] < csize[rb]) { tmp <- ra; ra <- rb; rb <- tmp }
    parent[rb] <<- ra
    csize[ra]  <<- csize[ra] + csize[rb]
    cid[ra]    <<- k
  }

  for (e in seq_len(nrow(el)))
    if (find(el[e, 1L]) != find(el[e, 2L])) link(el[e, 1L], el[e, 2L], ew[e])

  # join any disconnected components at height M
  roots <- unique(vapply(seq_len(n), find, integer(1L)))
  for (r in roots[-1L]) link(roots[1L], r, M)

  # iterative depth-first leaf order (avoids stack overflow on caterpillar trees)
  m  <- nrow(merge)
  res <- integer(n); ri <- 0L
  st  <- integer(2L * m + 2L); sp <- 1L; st[1L] <- m
  while (sp >= 1L) {
    nd <- st[sp]; sp <- sp - 1L
    if (nd < 0L) { ri <- ri + 1L; res[ri] <- -nd }
    else { st[sp + 1L] <- merge[nd, 2L]; st[sp + 2L] <- merge[nd, 1L]; sp <- sp + 2L }
  }

  labs <- if (!is.null(colnames(S))) colnames(S) else as.character(seq_len(n))
  hc   <- structure(
    list(merge = merge, height = height, order = res, labels = labs,
         method = "single", dist.method = "M - similarity"),
    class = "hclust"
  )
  ape::as.phylo(hc)
}


#' GetFaithPD (similarity-based)
#'
#' Calculates Faith's phylogenetic diversity using a pairwise feature
#' \strong{similarity} matrix from \code{GetSimMat()}. Converts to distance
#' internally (\code{1 - sim}) before clustering, which requires materializing
#' a full dense matrix; above 10,000 features it switches to a minimum spanning
#' tree, which stays sparse but uses single-linkage rather than average-linkage.
#'
#' @param feature Feature table with columns: id, feature, then sample columns
#' @param metadata Metadata table with sample and group columns
#' @param sim_matrix Feature similarity matrix (dense or sparse dgCMatrix)
#' @param threshold Numeric; detection threshold for presence (default: 0)
#' @param use_mst Logical; if TRUE, always use MST/single-linkage tree (works
#'   on sparse matrices of any size). If FALSE (default), average-linkage is
#'   used for n <= 10,000 and MST is used automatically (with a warning) for
#'   larger matrices.
#' @param use_fastcluster Logical; if TRUE, build the average-linkage tree with
#'   fastcluster::hclust() instead of stats::hclust(). Both are valid UPGMA
#'   implementations, but they break ties between equal distances differently.
#'   Chemical similarity scores are usually reported to three decimals, so ties
#'   are common and the two engines can return different trees. The default
#'   FALSE keeps results identical to earlier releases and to machines that do
#'   not have fastcluster installed.
#' @export
#' @return A data frame with sample, group, PD, SR, and value columns.
GetFaithPD <- function(feature, metadata, sim_matrix, threshold = 0, use_mst = FALSE,
                       use_fastcluster = FALSE) {
  .require_pkg("picante")
  .require_pkg("ape")
  n <- nrow(sim_matrix)

  # ------------------------------------------------------------------
  # Tree building: average-linkage for small n, MST for large n or on request
  # ------------------------------------------------------------------
  if (!use_mst && n <= 10000L) {
    dist_mat <- as.dist(1 - as.matrix(sim_matrix))
    hc <- if (isTRUE(use_fastcluster)) {
      .require_pkg("fastcluster")
      fastcluster::hclust(dist_mat, method = "average")
    } else {
      stats::hclust(dist_mat, method = "average")
    }
    tree <- ape::as.phylo(hc)
  } else {
    if (n > 10000L && !use_mst) {
      warning(
        "GetFaithPD: n = ", n, " features is too large for average-linkage clustering.\n",
        "Automatically switching to minimum spanning tree (single-linkage). ",
        "Results are valid but may differ from smaller datasets.\n",
        "To avoid this: use filter_mmo() to reduce feature count below 10,000.",
        call. = FALSE
      )
    } else {
      message("GetFaithPD: using minimum spanning tree (single-linkage) tree (use_mst = TRUE).")
    }
    tree <- sparse_single_phylo(sim_matrix, M = 1)
  }

  ids <- feature$id
  x <- as.matrix(feature[, -(1:2), drop = FALSE])
  # Faith PD uses presence/absence; threshold defines "present"
  x <- ifelse(!is.na(x) & x > threshold, 1, 0)
  feature_t <- t(x)
  colnames(feature_t) <- ids
  pd_result <- picante::pd(feature_t, tree)
  pd_result$sample <- rownames(pd_result)
  pd_result$group <- metadata$group[match(pd_result$sample, metadata$sample)]
  pd_result$value <- pd_result$PD
  rownames(pd_result) <- NULL
  return(pd_result)
}




#' GetAlphaDiversity
#'
#' Calculate alpha diversity for a \code{mmo} with flexible output modes.
#' Supported diversity modes:
#' - 'weighted'   : functional Hill number (GetFunctionalHillNumber)
#' - 'unweighted' : Hill numbers on abundances (GetHillNumbers)
#' - 'faith'      : Faith's phylogenetic diversity (GetFaithPD)
#' - 'richness'   : simple feature richness (GetRichness)
#'
#' Output modes control how samples are handled:
#' 1) 'sample_level'     : alpha per sample -> returns sample, group, value
#' 2) 'group_average'    : mean alpha per group (summarize sample_level) -> group, mean, sd, se, n, lwr, upr
#' 3) 'group_cumulative' : pooled gamma per group (pool samples within group) -> group, value
#' 4) 'rarefied_sample'  : sample-based rarefaction within group (subsample N samples, pool, compute) -> group, n_samples, mean, lwr, upr, n_perm
#'
#' NOTE: For outputs 3 and 4, pooling is performed by summing feature intensities across samples.
#'
#' @param mmo The \code{mmo} containing feature data and metadata
#' @param q The Hill number order controlling abundance sensitivity (default: 1).
#'   Only applies to \code{mode = "weighted"} and \code{mode = "unweighted"};
#'   ignored for \code{"richness"} and \code{"faith"}.
#'   \itemize{
#'     \item \code{q = 0} -- richness: all detected features count equally regardless of abundance.
#'     \item \code{q = 1} -- Shannon-type: features weighted proportionally to their relative abundance.
#'     \item \code{q = 2} -- Simpson-type: dominant (high-abundance) features weighted more strongly.
#'   }
#' @param normalization Abundance table to use. Options: 'None', 'Log', 'Meancentered', 'Z', 'PA'
#'   (default: 'None'). Using \code{'PA'} forces presence/absence regardless of \code{mode} or \code{q},
#'   effectively making every detected feature equally abundant before the Hill calculation.
#' @param mode The diversity metric to calculate. One of 'weighted', 'unweighted', 'faith', 'richness'
#'   (default: 'richness'). Use \code{q} to control abundance sensitivity for 'weighted' and 'unweighted'.
#' @param distance Feature dissimilarity metric: 'dreams', 'm2ds', or 'cosine' (default: 'dreams').
#'   Required for \code{mode = "weighted"} and \code{mode = "faith"}; ignored otherwise.
#' @param threshold Numeric threshold used to define metabolite presence (default: 0)
#' @param filter_id A boolean indicating whether to filter the feature data by a specific list (default: FALSE)
#' @param id_list A list of feature names to filter the feature data by, if filter_id is TRUE (default: NULL)
#' @param filter_group A boolean indicating whether to filter the feature data by a specific group list (default: FALSE)
#' @param group_list A list of groups to filter the feature data by, if filter_group is TRUE (default: NULL)
#'
#' @param output Output mode: 'sample_level', 'group_average', 'group_cumulative', or 'rarefied_sample'
#' @param group_col Column in mmo$metadata that defines groups (default: 'group')
#' @param sample_col Column in mmo$metadata that defines sample IDs (default: 'sample')
#' @param pool_method How to pool abundances when combining samples: 'sum' or 'mean' (default: 'sum')
#'
#' @param n_perm Integer; maximum number of permutations per rarefaction level (default: 200)
#' @param ci Numeric; confidence level (default: 0.95)
#' @param seed Optional integer seed for reproducibility (default: NULL)
#'
#' @return For output != 'rarefied_sample': a data.frame.
#'   For output = 'rarefied_sample': a list with:
#'   - summary: group-level rarefaction summary (mean, lwr, upr, n_perm_eff)
#'   - raw: permutation-level values for each group and n_samples
#' @export
# DEPRECATED: superseded by GetAlphaDiversity().
# GetAlphaDiversity_derep() uses GetDistanceMat() (.dissim slots, dissimilarity math).
# GetAlphaDiversity() uses GetSimMat() (.sim slots, similarity math).
GetAlphaDiversity_derep <- function(
    mmo,
    q = 1,
    normalization = "None",
    mode = "richness",
    distance = "dreams",
    threshold = 0,
    filter_id = FALSE,
    id_list = NULL,
    filter_group = FALSE,
    group_list = NULL,
    output = c("sample_level", "group_average", "group_cumulative", "rarefied_sample"),
    group_col = "group",
    sample_col = "sample",
    pool_method = c("sum", "mean"),
    n_perm = 200,
    ci = 0.95,
    seed = NULL
) {
  output <- match.arg(output)
  pool_method <- match.arg(pool_method)

  # -----------------------------
  # 0) Optional filtering
  # -----------------------------
  if (filter_id || filter_group) {
    mmo <- filter_mmo(mmo, id_list = id_list, group_list = group_list,
                      sample_col = sample_col, group_col = group_col)
  }

  # -----------------------------
  # 1) Metric engine (returns per-sample)
  # -----------------------------
  feature_all <- GetNormFeature(mmo, normalization = normalization)
  metadata_all <- mmo$metadata

  sample_cols_all <- colnames(feature_all)[-c(1, 2)]

  distance_matrix_all <- NULL
  if (mode %in% c("weighted", "faith")) {
    distance_matrix_all <- GetDistanceMat(mmo, distance = distance)
    keep_ids <- intersect(feature_all$id, rownames(distance_matrix_all))
    if (length(keep_ids) < 2) {
      stop("Need at least 2 overlapping features between feature data and distance matrix.")
    }
    feature_all <- feature_all[feature_all$id %in% keep_ids, , drop = FALSE]
    distance_matrix_all <- distance_matrix_all[keep_ids, keep_ids, drop = FALSE]
  }

  metric_engine <- function(feature_local, metadata_local, distance_local = NULL) {
    if (mode == "weighted") {
      if (is.null(distance_local)) {
        stop("Distance matrix required for weighted mode.")
      }
      GetFunctionalHillNumber_derep(
        feature = feature_local,
        metadata = metadata_local,
        distance_matrix = distance_local,
        q = q,
        threshold = threshold
      )
    } else if (mode == "unweighted") {
      GetHillNumbers(
        feature = feature_local,
        metadata = metadata_local,
        q = q,
        threshold = threshold
      )
    } else if (mode == "richness") {
      GetRichness(
        feature_data = feature_local,
        metadata = metadata_local,
        threshold = threshold
      )
    } else if (mode == "faith") {
      if (is.null(distance_local)) {
        stop("Distance matrix required for faith mode.")
      }
      GetFaithPD_derep(
        feature = feature_local,
        metadata = metadata_local,
        distance_matrix = distance_local,
        threshold = threshold
      )
    } else {
      stop("mode should be 'weighted', 'unweighted', 'richness', or 'faith'")
    }
  }

  # -----------------------------
  # 2) Standardize output to (sample, group, value)
  # -----------------------------
  to_sample_group_value <- function(df, mmo_ref) {
    if (!is.data.frame(df) || nrow(df) < 1) stop("Metric function returned empty output.")
    if (!("sample" %in% names(df))) {
      stop("Metric output must include a 'sample' column.")
    }
    if (!("group" %in% names(df))) {
      # map group from metadata
      md <- mmo_ref$metadata
      if (is.null(md) || !(sample_col %in% names(md)) || !(group_col %in% names(md))) {
        stop("Metric output lacks 'group' and cannot map it from mmo$metadata.")
      }
      df$group <- md[[group_col]][match(df$sample, md[[sample_col]])]
    }

    if ("richness" %in% names(df)) {
      value <- df$richness
    } else if ("value" %in% names(df)) {
      value <- df$value
    } else if ("alpha" %in% names(df)) {
      value <- df$alpha
    } else {
      stop("Metric output must include a numeric column named 'value' (or 'alpha'/'richness').")
    }

    data.frame(
      sample = as.character(df$sample),
      group  = as.character(df$group),
      value  = as.numeric(value),
      stringsAsFactors = FALSE
    )
  }

  # -----------------------------
  # 3) Pool feature columns helper
  # -----------------------------
  pool_feature_columns <- function(
      feature_in,
      metadata_in,
      samples = NULL,
      by_group = TRUE,
      pool_method = "sum"
  ) {
    sample_cols <- colnames(feature_in)[-c(1, 2)]
    x <- as.matrix(feature_in[, sample_cols, drop = FALSE])

    if (by_group) {
      groups <- unique(metadata_in[[group_col]])
      pooled <- matrix(NA_real_, nrow = nrow(x), ncol = length(groups))
      colnames(pooled) <- groups

      for (i in seq_along(groups)) {
        g <- groups[i]
        cols_g <- metadata_in[[sample_col]][metadata_in[[group_col]] == g]
        if (pool_method == "sum") {
          pooled[, i] <- rowSums(x[, cols_g, drop = FALSE], na.rm = TRUE)
        } else {
          pooled[, i] <- rowMeans(x[, cols_g, drop = FALSE], na.rm = TRUE)
        }
      }

      pooled_feature <- data.frame(
        feature_in[, c("id", "feature"), drop = FALSE],
        as.data.frame(pooled, check.names = FALSE),
        check.names = FALSE
      )
      pooled_metadata <- data.frame(sample = groups, group = groups, stringsAsFactors = FALSE)
    } else {
      if (is.null(samples) || length(samples) < 1) {
        stop("samples must be provided when by_group = FALSE.")
      }
      vals <- if (pool_method == "sum") {
        rowSums(x[, samples, drop = FALSE], na.rm = TRUE)
      } else {
        rowMeans(x[, samples, drop = FALSE], na.rm = TRUE)
      }
      pooled_feature <- data.frame(
        feature_in[, c("id", "feature"), drop = FALSE],
        pool = vals,
        check.names = FALSE
      )
      pooled_metadata <- data.frame(sample = "pool", group = "pool", stringsAsFactors = FALSE)
    }

    list(feature = pooled_feature, metadata = pooled_metadata)
  }

  # -----------------------------
  # Output 1) sample_level
  # -----------------------------
  if (output == "sample_level") {
    df <- metric_engine(feature_all, metadata_all, distance_matrix_all)
    return(to_sample_group_value(df, mmo))
  }

  # -----------------------------
  # Output 2) group_average
  # -----------------------------
  if (output == "group_average") {
    df <- metric_engine(feature_all, metadata_all, distance_matrix_all)
    sgv <- to_sample_group_value(df, mmo)

    alpha <- (1 - ci) / 2
    probs <- c(alpha, 1 - alpha)

    out <- sgv |>
      dplyr::group_by(.data$group) |>
      dplyr::summarise(
        mean = mean(.data$value, na.rm = TRUE),
        sd   = stats::sd(.data$value, na.rm = TRUE),
        n    = sum(!is.na(.data$value)),
        se   = .data$sd / sqrt(.data$n),
        lwr  = stats::quantile(.data$value, probs = probs[1], na.rm = TRUE),
        upr  = stats::quantile(.data$value, probs = probs[2], na.rm = TRUE),
        .groups = "drop"
      )
    return(as.data.frame(out))
  }

  # -----------------------------
  # Output 3) group_cumulative (pool within group, then compute metric)
  # -----------------------------
  if (output == "group_cumulative") {
    pooled_obj <- pool_feature_columns(
      feature_in = feature_all,
      metadata_in = metadata_all,
      by_group = TRUE,
      pool_method = pool_method
    )
    df <- metric_engine(pooled_obj$feature, pooled_obj$metadata, distance_matrix_all)
    sgv <- to_sample_group_value(df, mmo)

    # here sample names are the group labels (because pooling renames sample columns to groups)
    out <- data.frame(
      group = sgv$sample,
      value = sgv$value,
      stringsAsFactors = FALSE
    )
    return(out)
  }

  # -----------------------------
  # Output 4) rarefied_sample (curve up to max_k; sample-based)
  # -----------------------------
  if (output == "rarefied_sample") {

    if (!is.null(seed)) set.seed(seed)

    md <- metadata_all
    fd <- feature_all

    sample_cols <- sample_cols_all
    if (length(sample_cols) < 2) {
      stop("Need at least 2 aligned samples for rarefaction.")
    }

    group_map <- as.character(md$group[match(sample_cols, md$sample)])
    group_map[is.na(group_map) | group_map == ""] <- "ungrouped"
    groups <- unique(group_map)
    group_sizes <- table(group_map)
    single_sample_groups <- names(group_sizes)[group_sizes == 1]
    if (length(single_sample_groups) > 0) {
      warning(
        "rarefied_sample: group(s) with only 1 sample detected: ",
        paste(single_sample_groups, collapse = ", "),
        ". These groups will return only n_samples = 1 with n_perm_eff = 1.",
        call. = FALSE
      )
    }

    alpha <- (1 - ci) / 2
    probs <- c(alpha, 1 - alpha)

    # helper: maximum unique subsets for a given level k
    max_unique_subsets <- function(ng, k) {
      out <- suppressWarnings(choose(ng, k))
      if (!is.finite(out)) return(Inf)
      as.numeric(out)
    }

    res <- lapply(groups, function(g) {

      samp_g <- sample_cols[group_map == g]
      ng <- length(samp_g)
      k_vals <- seq_len(ng)
      mean_vals <- numeric(length(k_vals))
      lwr_vals <- numeric(length(k_vals))
      upr_vals <- numeric(length(k_vals))
      n_eff_vals <- integer(length(k_vals))
      raw_rows <- vector("list", length(k_vals))

      for (k in k_vals) {
        n_perm_eff_k <- min(n_perm, max_unique_subsets(ng, k))
        n_eff_vals[k] <- as.integer(n_perm_eff_k)

        vals_k <- numeric(n_eff_vals[k])
        for (p in seq_len(n_eff_vals[k])) {
          pick <- sample(samp_g, size = k, replace = FALSE)
          pooled_obj <- pool_feature_columns(
            feature_in = fd,
            metadata_in = md,
            samples = pick,
            by_group = FALSE,
            pool_method = pool_method
          )
          df_pool <- metric_engine(feature_local = pooled_obj$feature, metadata_local = pooled_obj$metadata, distance_local = distance_matrix_all)
          sgv_pool <- to_sample_group_value(df_pool, mmo)
          vals_k[p] <- sgv_pool$value[1]
        }

        raw_rows[[k]] <- data.frame(
          group = g,
          n_samples = k,
          perm_id = seq_len(n_eff_vals[k]),
          value = vals_k,
          stringsAsFactors = FALSE
        )

        mean_vals[k] <- mean(vals_k, na.rm = TRUE)
        lwr_vals[k]  <- stats::quantile(vals_k, probs = probs[1], na.rm = TRUE)
        upr_vals[k]  <- stats::quantile(vals_k, probs = probs[2], na.rm = TRUE)
      }

      list(
        summary = data.frame(
          group = g,
          n_samples = k_vals,
          mean = mean_vals,
          lwr  = lwr_vals,
          upr  = upr_vals,
          n_perm_eff = n_eff_vals,
          stringsAsFactors = FALSE
        ),
        raw = do.call(rbind, raw_rows)
      )
    })

    out_summary <- do.call(rbind, lapply(res, function(x) x$summary))
    out_raw <- do.call(rbind, lapply(res, function(x) x$raw))
    rownames(out_summary) <- NULL
    rownames(out_raw) <- NULL
    return(list(summary = out_summary, raw = out_raw))
  }

  stop("Unhandled output mode.")
}


#' GetAlphaDiversity (similarity-based)
#'
#' Identical to \code{GetAlphaDiversity_derep()} but uses \code{GetSimMat()}
#' to retrieve similarity matrices from \code{.sim} slots (created by
#' \code{AddChemSim()}) and passes them to the similarity-aware inner
#' functions (\code{GetFunctionalHillNumber()}, \code{GetFaithPD()}).
#'
#' @inheritParams GetAlphaDiversity_derep
#' @export
GetAlphaDiversity <- function(
    mmo,
    q = 1,
    normalization = "None",
    mode = "richness",
    distance = "dreams",
    threshold = 0,
    filter_id = FALSE,
    id_list = NULL,
    filter_group = FALSE,
    group_list = NULL,
    output = c("sample_level", "group_average", "group_cumulative", "rarefied_sample"),
    group_col = "group",
    sample_col = "sample",
    pool_method = c("sum", "mean"),
    n_perm = 200,
    ci = 0.95,
    seed = NULL,
    use_mst = FALSE,
    use_fastcluster = FALSE
) {
  output <- match.arg(output)
  pool_method <- match.arg(pool_method)

  # -----------------------------
  # 0) Optional filtering
  # -----------------------------
  if (filter_id || filter_group) {
    mmo <- filter_mmo(mmo, id_list = id_list, group_list = group_list,
                      sample_col = sample_col, group_col = group_col)
  }

  # -----------------------------
  # 1) Metric engine (returns per-sample)
  # -----------------------------
  feature_all <- GetNormFeature(mmo, normalization = normalization)
  metadata_all <- mmo$metadata

  sample_cols_all <- colnames(feature_all)[-c(1, 2)]

  distance_matrix_all <- NULL
  if (mode %in% c("weighted", "faith")) {
    distance_matrix_all <- GetSimMat(mmo, distance = distance)  # similarity .sim slots
    keep_ids <- intersect(feature_all$id, rownames(distance_matrix_all))
    if (length(keep_ids) < 2) {
      stop("Need at least 2 overlapping features between feature data and distance matrix.")
    }
    feature_all <- feature_all[feature_all$id %in% keep_ids, , drop = FALSE]
    distance_matrix_all <- distance_matrix_all[keep_ids, keep_ids, drop = FALSE]
  }

  metric_engine <- function(feature_local, metadata_local, distance_local = NULL) {
    if (mode == "weighted") {
      if (is.null(distance_local)) {
        stop("Distance matrix required for weighted mode.")
      }
      GetFunctionalHillNumber(
        feature = feature_local,
        metadata = metadata_local,
        sim_matrix = distance_local,
        q = q,
        threshold = threshold
      )
    } else if (mode == "unweighted") {
      GetHillNumbers(
        feature = feature_local,
        metadata = metadata_local,
        q = q,
        threshold = threshold
      )
    } else if (mode == "richness") {
      GetRichness(
        feature_data = feature_local,
        metadata = metadata_local,
        threshold = threshold
      )
    } else if (mode == "faith") {
      if (is.null(distance_local)) {
        stop("Distance matrix required for faith mode.")
      }
      GetFaithPD(
        feature = feature_local,
        metadata = metadata_local,
        sim_matrix = distance_local,
        threshold = threshold,
        use_mst = use_mst,
        use_fastcluster = use_fastcluster
      )
    } else {
      stop("mode should be 'weighted', 'unweighted', 'richness', or 'faith'")
    }
  }

  # -----------------------------
  # 2) Standardize output to (sample, group, value)
  # -----------------------------
  to_sample_group_value <- function(df, mmo_ref) {
    if (!is.data.frame(df) || nrow(df) < 1) stop("Metric function returned empty output.")
    if (!("sample" %in% names(df))) {
      stop("Metric output must include a 'sample' column.")
    }
    if (!("group" %in% names(df))) {
      md <- mmo_ref$metadata
      if (is.null(md) || !(sample_col %in% names(md)) || !(group_col %in% names(md))) {
        stop("Metric output lacks 'group' and cannot map it from mmo$metadata.")
      }
      df$group <- md[[group_col]][match(df$sample, md[[sample_col]])]
    }

    if ("richness" %in% names(df)) {
      value <- df$richness
    } else if ("value" %in% names(df)) {
      value <- df$value
    } else if ("alpha" %in% names(df)) {
      value <- df$alpha
    } else {
      stop("Metric output must include a numeric column named 'value' (or 'alpha'/'richness').")
    }

    data.frame(
      sample = as.character(df$sample),
      group  = as.character(df$group),
      value  = as.numeric(value),
      stringsAsFactors = FALSE
    )
  }

  # -----------------------------
  # 3) Pool feature columns helper
  # -----------------------------
  pool_feature_columns <- function(
      feature_in,
      metadata_in,
      samples = NULL,
      by_group = TRUE,
      pool_method = "sum"
  ) {
    sample_cols <- colnames(feature_in)[-c(1, 2)]
    x <- as.matrix(feature_in[, sample_cols, drop = FALSE])

    if (by_group) {
      groups <- unique(metadata_in[[group_col]])
      pooled <- matrix(NA_real_, nrow = nrow(x), ncol = length(groups))
      colnames(pooled) <- groups

      for (i in seq_along(groups)) {
        g <- groups[i]
        cols_g <- metadata_in[[sample_col]][metadata_in[[group_col]] == g]
        if (pool_method == "sum") {
          pooled[, i] <- rowSums(x[, cols_g, drop = FALSE], na.rm = TRUE)
        } else {
          pooled[, i] <- rowMeans(x[, cols_g, drop = FALSE], na.rm = TRUE)
        }
      }

      pooled_feature <- data.frame(
        feature_in[, c("id", "feature"), drop = FALSE],
        as.data.frame(pooled, check.names = FALSE),
        check.names = FALSE
      )
      pooled_metadata <- data.frame(sample = groups, group = groups, stringsAsFactors = FALSE)
    } else {
      if (is.null(samples) || length(samples) < 1) {
        stop("samples must be provided when by_group = FALSE.")
      }
      vals <- if (pool_method == "sum") {
        rowSums(x[, samples, drop = FALSE], na.rm = TRUE)
      } else {
        rowMeans(x[, samples, drop = FALSE], na.rm = TRUE)
      }
      pooled_feature <- data.frame(
        feature_in[, c("id", "feature"), drop = FALSE],
        pool = vals,
        check.names = FALSE
      )
      pooled_metadata <- data.frame(sample = "pool", group = "pool", stringsAsFactors = FALSE)
    }

    list(feature = pooled_feature, metadata = pooled_metadata)
  }

  # -----------------------------
  # Output 1) sample_level
  # -----------------------------
  if (output == "sample_level") {
    df <- metric_engine(feature_all, metadata_all, distance_matrix_all)
    return(to_sample_group_value(df, mmo))
  }

  # -----------------------------
  # Output 2) group_average
  # -----------------------------
  if (output == "group_average") {
    df <- metric_engine(feature_all, metadata_all, distance_matrix_all)
    sgv <- to_sample_group_value(df, mmo)

    alpha <- (1 - ci) / 2
    probs <- c(alpha, 1 - alpha)

    out <- sgv |>
      dplyr::group_by(.data$group) |>
      dplyr::summarise(
        mean = mean(.data$value, na.rm = TRUE),
        sd   = stats::sd(.data$value, na.rm = TRUE),
        n    = sum(!is.na(.data$value)),
        se   = .data$sd / sqrt(.data$n),
        lwr  = stats::quantile(.data$value, probs = probs[1], na.rm = TRUE),
        upr  = stats::quantile(.data$value, probs = probs[2], na.rm = TRUE),
        .groups = "drop"
      )
    return(as.data.frame(out))
  }

  # -----------------------------
  # Output 3) group_cumulative (pool within group, then compute metric)
  # -----------------------------
  if (output == "group_cumulative") {
    pooled_obj <- pool_feature_columns(
      feature_in = feature_all,
      metadata_in = metadata_all,
      by_group = TRUE,
      pool_method = pool_method
    )
    df <- metric_engine(pooled_obj$feature, pooled_obj$metadata, distance_matrix_all)
    sgv <- to_sample_group_value(df, mmo)

    # here sample names are the group labels (because pooling renames sample columns to groups)
    out <- data.frame(
      group = sgv$sample,
      value = sgv$value,
      stringsAsFactors = FALSE
    )
    return(out)
  }

  # -----------------------------
  # Output 4) rarefied_sample (curve up to max_k; sample-based)
  # -----------------------------
  if (output == "rarefied_sample") {

    if (!is.null(seed)) set.seed(seed)

    md <- metadata_all
    fd <- feature_all

    sample_cols <- sample_cols_all
    if (length(sample_cols) < 2) {
      stop("Need at least 2 aligned samples for rarefaction.")
    }

    group_map <- as.character(md$group[match(sample_cols, md$sample)])
    group_map[is.na(group_map) | group_map == ""] <- "ungrouped"
    groups <- unique(group_map)
    group_sizes <- table(group_map)
    single_sample_groups <- names(group_sizes)[group_sizes == 1]
    if (length(single_sample_groups) > 0) {
      warning(
        "rarefied_sample: group(s) with only 1 sample detected: ",
        paste(single_sample_groups, collapse = ", "),
        ". These groups will return only n_samples = 1 with n_perm_eff = 1.",
        call. = FALSE
      )
    }

    alpha <- (1 - ci) / 2
    probs <- c(alpha, 1 - alpha)

    # helper: maximum unique subsets for a given level k
    max_unique_subsets <- function(ng, k) {
      out <- suppressWarnings(choose(ng, k))
      if (!is.finite(out)) return(Inf)
      as.numeric(out)
    }

    res <- lapply(groups, function(g) {

      samp_g <- sample_cols[group_map == g]
      ng <- length(samp_g)
      k_vals <- seq_len(ng)
      mean_vals <- numeric(length(k_vals))
      lwr_vals <- numeric(length(k_vals))
      upr_vals <- numeric(length(k_vals))
      n_eff_vals <- integer(length(k_vals))
      raw_rows <- vector("list", length(k_vals))

      for (k in k_vals) {
        n_perm_eff_k <- min(n_perm, max_unique_subsets(ng, k))
        n_eff_vals[k] <- as.integer(n_perm_eff_k)

        vals_k <- numeric(n_eff_vals[k])
        for (p in seq_len(n_eff_vals[k])) {
          pick <- sample(samp_g, size = k, replace = FALSE)
          pooled_obj <- pool_feature_columns(
            feature_in = fd,
            metadata_in = md,
            samples = pick,
            by_group = FALSE,
            pool_method = pool_method
          )
          df_pool <- metric_engine(feature_local = pooled_obj$feature,
                                   metadata_local = pooled_obj$metadata,
                                   distance_local = distance_matrix_all)
          sgv_pool <- to_sample_group_value(df_pool, mmo)
          vals_k[p] <- sgv_pool$value[1]
        }

        raw_rows[[k]] <- data.frame(
          group = g,
          n_samples = k,
          perm_id = seq_len(n_eff_vals[k]),
          value = vals_k,
          stringsAsFactors = FALSE
        )

        mean_vals[k] <- mean(vals_k, na.rm = TRUE)
        lwr_vals[k]  <- stats::quantile(vals_k, probs = probs[1], na.rm = TRUE)
        upr_vals[k]  <- stats::quantile(vals_k, probs = probs[2], na.rm = TRUE)
      }

      list(
        summary = data.frame(
          group = g,
          n_samples = k_vals,
          mean = mean_vals,
          lwr  = lwr_vals,
          upr  = upr_vals,
          n_perm_eff = n_eff_vals,
          stringsAsFactors = FALSE
        ),
        raw = do.call(rbind, raw_rows)
      )
    })

    out_summary <- do.call(rbind, lapply(res, function(x) x$summary))
    out_raw <- do.call(rbind, lapply(res, function(x) x$raw))
    rownames(out_summary) <- NULL
    rownames(out_raw) <- NULL
    return(list(summary = out_summary, raw = out_raw))
  }

  stop("Unhandled output mode.")
}

#' RarefactionAUC
#'
#' This function calculates the rarefaction AUC for a given rarefied richness table
#'
#' @param rarefied_richness The rarefied richness object containing feature data, annotations, and pairwise comparisons
#' @param n_boot The number of bootstrap samples to use (default: 1000)
#' @param seed The seed to use for the bootstrap samples (default: 513)
#' @export
#' @return A list containing the rarefaction AUC for each group in the metadata, with columns for group and AUC
#' @examplesIf FALSE
#' RarefactionAUC(rarefied_richness, n_boot = 1000)
RarefactionAUC <- function(rarefied_richness, n_boot = 1000, seed = 513){
  if(!is.null(seed)) set.seed(seed)
  raw <- rarefied_richness$raw
  groups <- unique(raw$group)
  # Set AUC function
  trapz_auc <- function(x, y){
    sum(diff(x) * (head(y,-1) + tail(y,-1)) / 2)
  }
  # Set list for AUC values and summary
  auc_boot_list <- list()
  summary_list <- list()
  # Bootstrapping for each group
  for(g in groups){
    df_g <- raw[raw$group == g, ]
    n_samples_vec <- sort(unique(df_g$n_samples))
    auc_boot <- numeric(n_boot)
    for(i in seq_len(n_boot)){
      sampled <- sapply(n_samples_vec, function(n){
        vals <- df_g$value[df_g$n_samples == n]
        sample(vals, 1)
      })
      auc_boot[i] <- trapz_auc(n_samples_vec, sampled)
    }
    auc_df <- data.frame(
      group = g,
      auc = auc_boot
    )
    auc_boot_list[[g]] <- auc_df
    summary_list[[g]] <- data.frame(
      group = g,
      AUC_mean = mean(auc_boot),
      AUC_lwr = unname(quantile(auc_boot,0.025)),
      AUC_upr = unname(quantile(auc_boot,0.975))
    )
  }
  auc_boot_df <- do.call(rbind, auc_boot_list)
  summary_df <- do.call(rbind, summary_list)
  rownames(summary_df) <- NULL
  rownames(auc_boot_df) <- NULL
  return(list(
    auc_boot = auc_boot_df,
    summary = summary_df
  ))
}



#' GetSpecializationIndex
#'
#' This function calculates the specialization index for a given \code{mmo}, normalization method, and optional filtering by groups and features.
#'
#' @param mmo The \code{mmo} containing feature data and metadata
#' @param normalization The normalization method to use for feature data.
#'        Options are 'None', 'Log', 'Meancentered', or 'Z' (default: 'None')
#' @param filter_group A boolean indicating whether to filter the feature data by a specific group list (default: FALSE)
#' @param group_list A list of groups to filter the feature data by, if filter_group is TRUE (default: NULL)
#' @param filter_id A boolean indicating whether to filter the feature data by a specific list (default: FALSE)
#' @param id_list A list of feature names to filter the feature data by, if filter_id is TRUE (default: NULL)
#' @export
#' @return A data frame containing the specialization index for each group in the metadata, with columns for group and specialization index.
#' @examplesIf FALSE
#' specialization_index <- GetSpecializationIndex(mmo,
#'                                                normalization = 'None',
#'                                                filter_group = FALSE)
#' specialization_index <- GetSpecializationIndex(mmo,
#'                                                normalization = 'Z',
#'                                                filter_group = TRUE,
#'                                                group_list = c('Control', 'Treatment1'),
#'                                                filter_id = TRUE)
GetSpecializationIndex <- function(mmo, normalization = 'None', filter_group = FALSE, group_list = NULL, filter_id = FALSE, id_list = NULL){
  if (filter_id||filter_group){
    mmo <- filter_mmo(mmo, id_list = id_list, group_list = group_list)
  }
  metadata <- mmo$metadata
  feature <- GetNormFeature(mmo, normalization)
  # Get frequency matrix
  Pij <- feature[, -(1:2)]
  rownames(Pij) <- feature$feature
  Pij <- t(Pij)
  Pij <- Pij / rowSums(Pij)
  Pij[is.na(Pij)] <- 0

  # Get shannon diversity index
  # Pij.diversity <- Pij * log(Pij, base = 2)
  # Pij.diversity[is.na(Pij.diversity)] <- 0
  # Pij.diversity <- -rowSums(Pij.diversity)

  # Get specialization index
  Pi <- colSums(Pij)/nrow(Pij) # average frequency of each metabolite

  Si <- t((t(Pij)/Pi) * log((t(Pij)/Pi), base = 2)) # specialized degree
  Si[is.na(Si)] = 0
  Si <- colSums(Si)/ncol(Pij) # for each feature get average specialization index
  Pij.specialization <- colSums(t(Pij)*Si) # for each sample get specialization index
  # Retrieve group for each sample from metadata
  groups <- sapply(colnames(feature)[-c(1,2)], function(sample) {
    metadata$group[metadata$sample == sample]
  })
  output <- data.frame(group = groups, specialization = Pij.specialization)
}





#' GetBetaDiversity
#'
#' Calculate beta diversity (sample-to-sample dissimilarity) for a \code{mmo}.
#' Four methods are supported, differing in how they handle feature abundance and
#' structural relationships among compounds:
#'
#' - 'bray'    : Bray-Curtis dissimilarity. Feature-overlap only; no compound
#'               distance matrix required. Abundance sensitivity is controlled
#'               via \code{normalization}: use \code{'None'} for raw abundances,
#'               \code{'Log'} to compress dynamic range, or \code{'PA'} to reduce
#'               to presence/absence.
#'
#' - 'jaccard' : Jaccard dissimilarity. Presence/absence only by definition.
#'               The \code{normalization} argument is ignored; PA data are always
#'               used internally. A warning is issued if a non-PA normalization
#'               is supplied.
#'
#' - 'CSCS'    : Chemical Structural and Compositional Similarity. Incorporates
#'               pairwise compound distances so that samples sharing structurally
#'               related (but not identical) features can still be considered
#'               similar. Requires \code{distance}. Abundance sensitivity is
#'               controlled via \code{normalization} (same options as bray).
#'
#' - 'Gen.Uni' : Generalized UniFrac. Partitions beta diversity across branches
#'               of a compound dendrogram. Requires \code{distance}. Returns a
#'               named list of three distance matrices, one per abundance-weighting
#'               level (\code{d_0}, \code{d_0.5}, \code{d_1}). Choose the one
#'               appropriate for your analysis (see Details).
#'
#' @details
#' \strong{Choosing an abundance-weighting level for Gen.Uni:}
#' \itemize{
#'   \item \code{d_0}   -- presence/absence only; equivalent to unweighted UniFrac.
#'                         Use when detection (not intensity) is the signal of interest.
#'   \item \code{d_0.5} -- balanced weighting. Moderates the
#'                         influence of highly abundant features without ignoring abundance.
#'   \item \code{d_1}   -- fully abundance-weighted. Dominant features drive the distance.
#'                         Use when peak intensity is a reliable biological signal.
#' }
#' Access a specific matrix with e.g. \code{result[["d_0.5"]]}.
#'
#' @param mmo The \code{mmo} containing feature data and metadata
#' @param method Beta diversity method: 'Gen.Uni', 'bray', 'jaccard', or 'CSCS' (default: 'Gen.Uni')
#' @param normalization Abundance table to use. Options: 'None', 'Log', 'Meancentered', 'Z', 'PA'
#'   (default: 'None'). Ignored for 'jaccard', which always uses PA.
#'   For 'bray' and 'CSCS', this is the primary lever for controlling abundance sensitivity.
#' @param distance Feature dissimilarity metric: 'dreams', 'm2ds', or 'cosine'.
#'   Required for 'Gen.Uni' and 'CSCS'; ignored for 'bray' and 'jaccard'.
#' @param filter_id A boolean indicating whether to filter the feature data by a specific list (default: FALSE)
#' @param id_list A list of feature names to filter the feature data by, if filter_id is TRUE (default: NULL)
#' @param filter_group A boolean indicating whether to filter the feature data by a specific group list (default: FALSE)
#' @param group_list A list of groups to filter the feature data by, if filter_group is TRUE (default: NULL)
#' @param scale_dissim Boolean; whether to scale the feature distance matrix to between 0,1 (default: TRUE)
#' @return For 'bray', 'jaccard', 'CSCS': a symmetric sample-by-sample distance matrix (matrix).
#'   For 'Gen.Uni': a named list with three distance matrices named \code{d_0} (presence/absence),
#'   \code{d_0.5} (balanced), and \code{d_1} (fully abundance-weighted).
#' @export
#' @examplesIf FALSE
#' beta_diversity <- GetBetaDiversity(mmo, method = 'Gen.Uni',
#'  normalization = 'None', distance = 'dreams', filter_id = FALSE)
#' # Access a specific weighting level:
#' d_balanced <- beta_diversity[["d_0.5"]]
#' beta_diversity <- GetBetaDiversity(mmo, method = 'bray',
#'  normalization = 'Z', filter_id = TRUE, id_list = Glucosinolates,
#'  filter_group = TRUE, group_list = c('Control', 'Treatment1'))
# DEPRECATED: superseded by GetBetaDiversity().
# GetBetaDiversity_derep() uses GetDistanceMat() (.dissim slots).
# GetBetaDiversity() uses GetSimMat() (.sim slots) with similarity-aware math:
#   CSCS uses S directly (no 1-D conversion); Gen.Uni converts S->D only
#   inside hclust (size-guarded for n > 10,000).
GetBetaDiversity_derep <- function(mmo, method = 'Gen.Uni', normalization = 'None', distance = NULL, filter_id = FALSE, id_list = NULL, filter_group = FALSE, group_list = NULL, scale_dissim = TRUE){
  if (filter_id||filter_group){
    mmo <- filter_mmo(mmo, id_list = id_list, group_list = group_list)
  }

  # Methods that require a compound distance matrix
  distance_methods <- c('Gen.Uni', 'CSCS')

  if (method %in% distance_methods) {
    if (is.null(distance)) stop(paste("method =", method, "requires a 'distance' argument."))
    dist_mat <- GetDistanceMat(mmo, distance = distance)
    scaled_dissimilarity <- if (scale_dissim) dist_mat / max(dist_mat) else dist_mat
  }

  # Build tree only for Gen.Uni
  if (method == 'Gen.Uni') {
    compound_tree <- ape::as.phylo(hclust(as.dist(scaled_dissimilarity), method = "average"))
  }

  # Get feature matrix
  metadata <- mmo$metadata
  feature <- GetNormFeature(mmo, normalization)

  if (method %in% distance_methods) {
    feature <- feature |> filter(.data$id %in% colnames(scaled_dissimilarity))
    relative_proportions <- apply(feature[, -(1:2)], 2, function(x) x / sum(x))
    rownames(relative_proportions) <- feature$id
    relative_proportions <- relative_proportions[rownames(scaled_dissimilarity), ]
    relative_proportions <- t(relative_proportions)
  }

  if (method == 'Gen.Uni') {
    guni <- GUniFrac::GUniFrac(relative_proportions, compound_tree, alpha = c(0, 0.5, 1), verbose = TRUE)
    beta_div <- list(
      d_0   = guni$unifracs[, , "d_0"],
      d_0.5 = guni$unifracs[, , "d_0.5"],
      d_1   = guni$unifracs[, , "d_1"]
    )
    message(
      "Gen.Uni returned three distance matrices with different abundance-weighting levels:\n",
      "  d_0   : presence/absence only (unweighted)\n",
      "  d_0.5 : balanced weighting\n",
      "  d_1   : fully abundance-weighted\n",
      "Access one with e.g. result[[\"d_0.5\"]] before passing to NMDSplot() or PCoAplot()."
    )
  } else if (method == 'bray') {
    relative_proportions <- apply(feature[, -(1:2)], 2, function(x) x / sum(x))
    beta_div <- as.matrix(vegan::vegdist(t(relative_proportions), method = 'bray'))
  } else if (method == 'jaccard') {
    if (normalization != 'PA') {
      warning(
        "Jaccard dissimilarity requires presence/absence data. ",
        "The 'normalization' argument ('", normalization, "') is ignored; ",
        "'PA' is always used. To silence this warning, pass normalization = 'PA'.",
        call. = FALSE
      )
    }
    feature <- GetNormFeature(mmo, normalization = 'PA')
    beta_div <- as.matrix(vegan::vegdist(t(feature[, -(1:2)]), method = 'jaccard'))
  } else if (method == 'CSCS') {
    CSS <- 1 - scaled_dissimilarity  # already computed above
    if(!all(diag(CSS) == 1)) print("Warning Diag not equal to 1")
    q.feature <- GetNormFeature(mmo, normalization = normalization) |> filter(.data$id %in% colnames(CSS))
    q.feature <- q.feature[match(colnames(CSS), q.feature$id), ]
    q.mat <- as.matrix(q.feature[, -(1:2), drop = FALSE])
    rownames(q.mat) <- q.feature$id
    relative_proportions <- sweep(q.mat, 2, colSums(q.mat), "/")
    if (!identical(rownames(relative_proportions), rownames(CSS))) {
      stop("Feature alignment failed: rownames(relative_proportions) != rownames(CSS)")
    }
    if (!identical(rownames(relative_proportions), colnames(CSS))) {
      stop("Feature alignment failed: rownames(relative_proportions) != colnames(CSS)")
    }
    CSCS_all <- t(relative_proportions) %*% CSS %*% relative_proportions
    sample_names <- colnames(relative_proportions)
    n_samples <- length(sample_names)
    CSCS_matrix <- matrix(NA, nrow = n_samples, ncol = n_samples)
    rownames(CSCS_matrix) <- sample_names
    colnames(CSCS_matrix) <- sample_names
    for (i in 1:n_samples) {
      for (j in 1:n_samples) {
        CSCS_matrix[i,j] <- CSCS_all[i, j] / max(CSCS_all[i, i], CSCS_all[j, j])
      }
    }
    beta_div <- 1 - CSCS_matrix
  } else {
    stop("Invalid method. Please use 'Gen.Uni', 'bray', 'jaccard', or 'CSCS'")
  }

  return(beta_div)
}


#' GetBetaDiversity (similarity-based)
#'
#' Identical to \code{GetBetaDiversity_derep()} but retrieves pairwise
#' \strong{similarity} matrices via \code{GetSimMat()} (.sim slots).
#' Key differences from the _derep version:
#' \itemize{
#'   \item CSCS: uses the similarity matrix directly as CSS (no \code{1 - D}
#'         conversion needed), because \code{CSS = S} when similarities are stored.
#'   \item Gen.Uni: converts \code{S -> D} internally (\code{1 - S}) only for
#'         the \code{hclust()} call; guarded by a size check (n > 10,000 errors).
#'   \item bray / jaccard: unchanged (no distance matrix used).
#' }
#'
#' @inheritParams GetBetaDiversity_derep
#' @inheritParams GetFaithPD
#' @export
GetBetaDiversity <- function(mmo, method = 'Gen.Uni', normalization = 'None', distance = NULL, filter_id = FALSE, id_list = NULL, filter_group = FALSE, group_list = NULL, scale_dissim = TRUE, use_mst = FALSE, use_fastcluster = FALSE){
  if (filter_id||filter_group){
    mmo <- filter_mmo(mmo, id_list = id_list, group_list = group_list)
  }

  # Methods that require a compound similarity matrix
  distance_methods <- c('Gen.Uni', 'CSCS')

  if (method %in% distance_methods) {
    if (is.null(distance)) stop(paste("method =", method, "requires a 'distance' argument."))
    sim_mat <- GetSimMat(mmo, distance = distance)
    scaled_similarity <- if (scale_dissim) sim_mat / max(sim_mat) else sim_mat
  }

  # Build tree only for Gen.Uni
  if (method == 'Gen.Uni') {
    n_feat <- nrow(sim_mat)
    if (!use_mst && n_feat <= 10000L) {
      dist_s <- as.dist(1 - as.matrix(scaled_similarity))
      hc_s <- if (isTRUE(use_fastcluster)) {
        .require_pkg("fastcluster")
        fastcluster::hclust(dist_s, method = "average")
      } else {
        stats::hclust(dist_s, method = "average")
      }
      compound_tree <- ape::as.phylo(hc_s)
    } else {
      if (n_feat > 10000L && !use_mst) {
        warning(
          "GetBetaDiversity(Gen.Uni): n = ", n_feat, " features is too large for ",
          "average-linkage clustering.\n",
          "Automatically switching to minimum spanning tree (single-linkage). ",
          "Results are valid but may differ from smaller datasets.\n",
          "To avoid this: use filter_mmo() to reduce feature count below 10,000.",
          call. = FALSE
        )
      } else {
        message("GetBetaDiversity(Gen.Uni): using minimum spanning tree (single-linkage) tree (use_mst = TRUE).")
      }
      compound_tree <- sparse_single_phylo(scaled_similarity, M = 1)
    }
  }

  # Get feature matrix
  metadata <- mmo$metadata
  feature <- GetNormFeature(mmo, normalization)

  if (method %in% distance_methods) {
    feature <- feature |> filter(.data$id %in% colnames(scaled_similarity))
    relative_proportions <- apply(feature[, -(1:2)], 2, function(x) x / sum(x))
    rownames(relative_proportions) <- feature$id
    relative_proportions <- relative_proportions[rownames(scaled_similarity), ]
    relative_proportions <- t(relative_proportions)
  }

  if (method == 'Gen.Uni') {
    guni <- GUniFrac::GUniFrac(relative_proportions, compound_tree, alpha = c(0, 0.5, 1), verbose = TRUE)
    beta_div <- list(
      d_0   = guni$unifracs[, , "d_0"],
      d_0.5 = guni$unifracs[, , "d_0.5"],
      d_1   = guni$unifracs[, , "d_1"]
    )
    message(
      "Gen.Uni returned three distance matrices with different abundance-weighting levels:\n",
      "  d_0   : presence/absence only (unweighted)\n",
      "  d_0.5 : balanced weighting\n",
      "  d_1   : fully abundance-weighted\n",
      "Access one with e.g. result[[\"d_0.5\"]] before passing to NMDSplot() or PCoAplot()."
    )
  } else if (method == 'bray') {
    relative_proportions <- apply(feature[, -(1:2)], 2, function(x) x / sum(x))
    beta_div <- as.matrix(vegan::vegdist(t(relative_proportions), method = 'bray'))
  } else if (method == 'jaccard') {
    if (normalization != 'PA') {
      warning(
        "Jaccard dissimilarity requires presence/absence data. ",
        "The 'normalization' argument ('", normalization, "') is ignored; ",
        "'PA' is always used. To silence this warning, pass normalization = 'PA'.",
        call. = FALSE
      )
    }
    feature <- GetNormFeature(mmo, normalization = 'PA')
    beta_div <- as.matrix(vegan::vegdist(t(feature[, -(1:2)]), method = 'jaccard'))
  } else if (method == 'CSCS') {
    # With similarity storage, CSS = scaled_similarity directly (no 1-D conversion)
    CSS <- scaled_similarity
    # Matrix::diag() reads the diagonal without expanding a sparse matrix to
    # dense, which as.matrix() would do purely to check n values.
    if (!all(Matrix::diag(CSS) == 1)) message("Warning: similarity diagonal is not all 1")
    q.feature <- GetNormFeature(mmo, normalization = normalization) |> filter(.data$id %in% colnames(CSS))
    q.feature <- q.feature[match(colnames(CSS), q.feature$id), ]
    q.mat <- as.matrix(q.feature[, -(1:2), drop = FALSE])
    rownames(q.mat) <- q.feature$id
    relative_proportions <- sweep(q.mat, 2, colSums(q.mat), "/")
    if (!identical(rownames(relative_proportions), rownames(CSS))) {
      stop("Feature alignment failed: rownames(relative_proportions) != rownames(CSS)")
    }
    if (!identical(rownames(relative_proportions), colnames(CSS))) {
      stop("Feature alignment failed: rownames(relative_proportions) != colnames(CSS)")
    }
    CSCS_all <- t(relative_proportions) %*% CSS %*% relative_proportions
    sample_names <- colnames(relative_proportions)
    n_samples <- length(sample_names)
    CSCS_matrix <- matrix(NA, nrow = n_samples, ncol = n_samples)
    rownames(CSCS_matrix) <- sample_names
    colnames(CSCS_matrix) <- sample_names
    for (i in 1:n_samples) {
      for (j in 1:n_samples) {
        CSCS_matrix[i,j] <- CSCS_all[i, j] / max(CSCS_all[i, i], CSCS_all[j, j])
      }
    }
    beta_div <- 1 - CSCS_matrix
  } else {
    stop("Invalid method. Please use 'Gen.Uni', 'bray', 'jaccard', or 'CSCS'")
  }

  return(beta_div)
}


#' NMDSplot
#'
#' This function generates a Non-metric Multidimensional Scaling (NMDS) plot based on the provided beta diversity distance matrix.
#' It also performs PERMANOVA analysis to assess the significance of group differences and saves the results to CSV files.
#' @param mmo The \code{mmo} containing metadata
#' @param betadiv The beta diversity distance matrix, output of GetBetaDiversity()
#' @param outdir The outdir for the output files
#' @param width The width of the output NMDS plot (default: 6)
#' @param height The height of the output NMDS plot (default: 6)
#' @param color A vector of colors for the groups in the plot
#' @param save_output A logical value indicating whether to save the output plot (default: TRUE)
#' @return A list containing the NMDS plot, the NMDS coordinates, and the PERMANOVA results
#' @export
#' @examplesIf FALSE
#' beta_diversity <- GetBetaDiversity(mmo, method = 'Gen.Uni',
#'  normalization = 'None', distance = 'dreams', filter_id = FALSE)
#' # Use method = 'bray' or 'jaccard' if you want to use just feature abundance
#' # without considering feature spectral dissimilarity
#' NMDSplot(mmo, betadiv = beta_diversity, outdir = 'output/NMDS', width = 6, height = 6)
NMDSplot <- function(mmo, betadiv, outdir, width = 6, height = 6, color = NULL, save_output = TRUE){
  .require_pkg("vegan")
  .require_pkg("ggrepel")
  metadata <- mmo$metadata
  nmds <- vegan::metaMDS(betadiv, k = 2, try = 50, trymax = 100)

  # If color is NULL, set default colors
  if (is.null(color)) {
    groups <- unique(metadata$group)
    colors <- setNames(grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(groups)), groups)
  }
  # Extract NMDS coordinates
  nmds_coords <- as.data.frame(vegan::scores(nmds, display = "sites"))
  groups <- c()
  for (row in rownames(nmds_coords)) {
    groups <- append(groups, metadata[metadata$sample == row, ]$group)
  }
  nmds_coords$group <- groups

  # Plot NMDS
  plot <- ggplot(nmds_coords, aes(x = .data$NMDS1, y = .data$NMDS2, color = .data$group)) +
    geom_point(size = 3) +
    #geom_text_repel(aes(label = group), size = 3) +
    theme_classic() +
    stat_ellipse(level = 0.90) +
    labs(x = "NMDS1", y = "NMDS2") +
    scale_color_manual(values = color) +
    theme(legend.position = "right")
  plot
  permanova <- permanova_stat(betadiv, mmo$metadata, mode = 'distance')
  if (save_output){
    ggsave(paste0(outdir, '_NMDS.pdf'), height = height, width = width)
    readr::write_csv(permanova$permanova_res, paste0(outdir, '_permanova_results.csv'))
    readr::write_csv(as.data.frame(permanova$pairwise_raw), paste0(outdir, '_pairwise_permanova_results.csv'))
    readr::write_csv(as.data.frame(permanova$pairwise_p_matrix), paste0(outdir, '_pairwise_permanova_pvalue_matrix.csv'))
    readr::write_csv(as.data.frame(permanova$pairwise_F_matrix), paste0(outdir, '_pairwise_permanova_F_matrix.csv'))
    readr::write_csv(as.data.frame(permanova$pairwise_R2_matrix), paste0(outdir, '_pairwise_permanova_R2_matrix.csv'))
  }
  return(list(plot = plot, df = nmds_coords, permanova = permanova))
}

#' PCoAplot
#'
#' This function generates a Principal Coordinates Analysis (PCoA) plot based on the provided beta diversity distance matrix.
#' It also performs PERMANOVA analysis to assess the significance of group differences and saves the
#' results to CSV files.
#' @param mmo The \code{mmo} containing metadata
#' @param betadiv The beta diversity distance matrix, output of GetBetaDiversity
#' @param outdir The prefix for the output files
#' @param width The width of the output PCoA plot (default: 6
#' @param height The height of the output PCoA plot (default: 6)
#' @param color A vector of colors for the points in the plot
#' @param save_output A logical value indicating whether to save the output plot (default: TRUE)
#' @return A list containing the PCoA plot, the PCoA coordinates, and the PERMANOVA results
#' @export
#' @examplesIf FALSE
#' beta_diversity <- GetBetaDiversity(mmo, method = 'Gen.Uni',
#'  normalization = 'None', distance = 'dreams', filter_id = FALSE)
#' PCoAplot(mmo, betadiv = beta_diversity, outdir = 'output/PCoA', width = 6, height = 6)
PCoAplot <- function(mmo, betadiv, outdir, width = 6, height = 6, color = NULL, save_output = TRUE){
  .require_pkg('ape')
  metadata <- mmo$metadata
  pcoa_res <- ape::pcoa(betadiv)
  pcoa_coords <- as.data.frame(pcoa_res$vectors[, 1:2])
  colnames(pcoa_coords) <- c("PCoA1", "PCoA2")
  pcoa_coords$group <- metadata$group[match(rownames(pcoa_coords), metadata$sample)]

  # If color is NULL, set default colors
  if (is.null(color)) {
    groups <- unique(metadata$group)
    colors <- setNames(grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(groups)), groups)
  }
  plot <- ggplot(pcoa_coords, aes(x = .data$PCoA1, y = .data$PCoA2, color = .data$group)) +
    geom_point(size = 3) +
    stat_ellipse(level = 0.90) +
    theme_classic() +
    labs(x = "PCoA1", y = "PCoA2") +
    scale_color_manual(values = color) +
    theme(legend.position = "right")
  plot
  permanova <- permanova_stat(betadiv, mmo$metadata, mode = 'distance')
  if (save_output){
    ggsave(paste0(outdir, '_PCoA.pdf'), height = height, width = width)
    readr::write_csv(permanova$permanova_res, paste0(outdir, '_permanova_results.csv'))
    readr::write_csv(as.data.frame(permanova$pairwise_raw), paste0(outdir, '_pairwise_permanova_results.csv'))
    readr::write_csv(as.data.frame(permanova$pairwise_p_matrix), paste0(outdir, '_pairwise_permanova_pvalue_matrix.csv'))
    readr::write_csv(as.data.frame(permanova$pairwise_F_matrix), paste0(outdir, '_pairwise_permanova_F_matrix.csv'))
    readr::write_csv(as.data.frame(permanova$pairwise_R2_matrix), paste0(outdir, '_pairwise_permanova_R2_matrix.csv'))
  }
  return(list(plot = plot, df = pcoa_coords, permanova = permanova))
}

#' CalculateGroupBetaDistance
#'
#' This function calculates the beta diversity distance between a reference group and other groups in the metadata.
#'
#' @param mmo The \code{mmo} containing feature data and metadata
#' @param beta_div The beta diversity distance matrix, output of GetBetaDiversity()
#' @param reference_group The name of the reference group to compare against
#' @param groups A vector of group names from the metadata to calculate distances for
#' @return A data frame containing the group names, sample names, and their corresponding beta diversity distances from the reference group.
#' @export
#' @examplesIf FALSE
#' beta_diversity <- GetBetaDiversity(mmo, method = 'Gen.Uni',
#'  normalization = 'None', distance = 'dreams', filter_id = FALSE)
#' group_distances <- CalculateGroupBetaDistance(mmo, beta_div = beta_diversity,
#'  reference_group = 'Control', groups = c('Control', 'Treatment1', 'Treatment2'))
CalculateGroupBetaDistance <- function(mmo, beta_div, reference_group, groups) {
  metadata <- mmo$metadata
  distances <- data.frame(group = character(), distance = numeric())

  for (group in groups) {
      group_samples <- metadata |> filter(group == !!group) |> pull(sample)
      reference_samples <- metadata |> filter(group == !!reference_group) |> pull(sample)

      for (sample in group_samples) {
        for (ref_sample in reference_samples) {
          distance <- beta_div[sample, ref_sample]
          distances <- rbind(distances, data.frame(group = group,sample = sample, distance = distance))
        }

    }
  }

  return(distances)
}

#' Save entire mmo object to a file (RDS)
#'
#' @param mmo The \code{mmo} to save
#' @param file File path to write (default: "mmo.rds")
#' @param compress Compression type passed to saveRDS ("gzip", "bzip2", "xz", or logical) (default: "xz")
#' @param include_session Logical; if TRUE attach sessionInfo() as an attribute to the saved object (default: TRUE)
#' @return Invisibly returns the file path
#' @export
SaveMMO <- function(mmo, file = "mmo.rds", compress = "xz", include_session = TRUE) {
  if (missing(mmo) || !is.list(mmo)) stop("mmo must be a list-like object")
  if (include_session) {
    attr(mmo, "saved_session_info") <- utils::sessionInfo()
  }
  saveRDS(mmo, file = file, compress = compress)
  message("Saved mmo to: ", file)
  invisible(file)
}


#' Load a \code{mmo} previously saved with SaveMMO
#'
#' This function returns the loaded \code{mmo} (visible return). By default it prints basic
#' information about the R version and recorded packages that were present when the object
#' was saved.
#'
#' @param file Path to an RDS file created with SaveMMO
#' @param check_session Logical; if TRUE and save-time session info is present, print a short summary (default: TRUE)
#' @param verbose Logical; print messages about saved session info when available (default: TRUE)
#' @return The loaded \code{mmo}
#' @export
LoadMMO <- function(file, check_session = TRUE, verbose = TRUE) {
  if (!file.exists(file)) stop("File not found: ", file)
  mmo <- readRDS(file)
  if (check_session && !is.null(attr(mmo, "saved_session_info"))) {
    saved_si <- attr(mmo, "saved_session_info")
    if (verbose) {
      message("mmo was saved with R: ", saved_si$R.version$version.string)
      pkgs <- names(saved_si$otherPkgs)
      if (length(pkgs) > 0) message("Top packages at save-time: ", paste(head(pkgs, 20), collapse = ", "))
    }
  }
  return(mmo)
}


#' Print method for \code{mmo}s
#' Provides a clean, human-readable overview of an `mmo` list object instead of
#' dumping the entire list when the object is printed in the console.
#' @param x An `mmo` object (a list with components such as `feature_data`,
#'   `metadata`, `pairwise`, etc.)
#' @param ... Additional arguments passed to other print methods (unused).
#' @return Invisibly returns `x` unchanged.
#' @export
print.mmo <- function(x, ...) {
  cat("MMO object\n")
  # Features
  if (!is.null(x$feature_data)) {
    cat("  Feature number: ",
        nrow(x$feature_data), "\n", sep = "")
  }

  # Samples & groups
  if (!is.null(x$metadata)) {
    n_samples <- nrow(x$metadata)
    n_groups  <- length(unique(x$metadata$group))
    cat("  ", n_samples, " samples in ", n_groups, " groups\n", sep = "")
  }

  # Components present
  cat("  MMO object contains: ",
      paste(names(x), collapse = ", "),
      "\n", sep = "")

  invisible(x)
}

#' HCplot
#'
#' Hierarchical clustering of samples from a precomputed beta-diversity matrix,
#' plotted as a phylogram with tip labels colored by \code{mmo$metadata$group}.
#'
#' @param mmo The \code{mmo}
#' @param betadiv Beta diversity distance matrix, output of GetBetaDiversity()
#' @param outdir Output file prefix (e.g., "output/HC"). PNG saved here if save_output=TRUE.
#' @param hclust_method hclust linkage method (default: "average"; alternatives: "complete","ward.D2")
#' @param color Named vector of colors for groups. If NULL, uses Set3 palette.
#' @param tip_label_size Tip label text size (default: 2.5)
#' @param width Plot width in inches (default: 10)
#' @param height Plot height in inches (default: 7)
#' @param save_output Whether to save the plot (default: TRUE)
#' @return A list containing: plot (ggplot object), hc (hclust), phy (phylo)
#' @export
#' @examplesIf FALSE
#' bet <- GetBetaDiversity(mmo, method = 'bray', normalization = 'None')
#' hc <- HCplot(mmo, betadiv = bet, outdir = "output/HC_bray")
#' hc$plot
HCplot <- function(
    mmo,
    betadiv,
    outdir,
    hclust_method = "average",
    color = NULL,
    tip_label_size = 2.5,
    width = 10,
    height = 7,
    save_output = TRUE
){
  .require_pkg("ape")
  .require_pkg("ggtree")

  metadata <- mmo$metadata

  if (is.null(color)) {
    groups <- unique(metadata$group)
    color <- setNames(grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(groups)), groups)
  }

  # Cluster
  d  <- as.dist(betadiv)
  hc <- stats::hclust(d, method = hclust_method)
  tr <- ape::as.phylo(hc)

  # Map group to tips
  tip_df <- data.frame(
    label = tr$tip.label,
    group = metadata$group[match(tr$tip.label, metadata$sample)],
    stringsAsFactors = FALSE
  )

  tip_depth <- max(ape::node.depth.edgelength(tr))

  plot <- ggtree::`%<+%`(ggtree::ggtree(tr, layout = "rectangular"), tip_df) +
    ggtree::geom_tippoint(ggplot2::aes(color = .data$group), size = 2) +
    ggtree::geom_tiplab(ggplot2::aes(color = .data$group),
                        size = tip_label_size, hjust = -0.05, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = color, name = "group") +
    ggplot2::labs(title = paste0("Hierarchical clustering (", hclust_method, ")")) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "right",
                   plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::xlim(NA, tip_depth * 1.5)

  if (save_output) {
    ggsave(paste0(outdir, '_HC.pdf'), plot = plot, width = width, height = height)
  }

  return(list(plot = plot, hc = hc, phy = tr))
}


#' FeatureDendrogram
#'
#' Build a feature-to-feature dendrogram from a dissimilarity matrix stored in
#' the \code{mmo}. Optionally collapses within-group distances for ion identity
#' or correlation groups so that adducts/isotopes of the same compound cluster
#' together before chemical distance drives the higher-level topology.
#'
#' @details
#' \strong{Ion identity options:}
#' \itemize{
#'   \item \code{"none"} -- use the raw dissimilarity matrix as-is.
#'   \item \code{"correlation"} -- features sharing the same MZmine correlation
#'     group (typically co-eluting features likely derived from the same
#'     compound) are assigned \code{within_group_dist} before clustering.
#'   \item \code{"ion_identity_network"} -- features linked by MZmine ion
#'     identity networking (confirmed adducts / isotopes) are assigned
#'     \code{within_group_dist} before clustering.
#' }
#' In both grouped modes the constraint is applied symmetrically and the
#' diagonal is forced to 0 after the assignment.
#'
#' @param mmo \code{mmo}. Must contain the requested dissimilarity matrix
#'   (added via \code{AddChemDist()} or \code{AddCustomDist()}) and, for
#'   ion identity modes, \code{mmo$feature_info}.
#' @param distance Name of the dissimilarity matrix to use (default: \code{"dreams"}).
#'   Passed to \code{GetDistanceMat()}; supports built-in and custom names.
#' @param features Optional character vector of feature IDs to include.
#'   \code{NULL} (default) uses all features in the distance matrix.
#' @param method hclust linkage method (default: \code{"average"}).
#'   Other sensible choices: \code{"complete"}, \code{"ward.D2"}.
#' @param ion_identity One of \code{"none"}, \code{"correlation"},
#'   \code{"ion_identity_network"} (default: \code{"none"}).
#' @param corr_col Column in \code{mmo$feature_info} containing MZmine
#'   correlation group IDs (default: \code{"feature_group"}).
#'   Used when \code{ion_identity = "correlation"}.
#' @param iin_col Column in \code{mmo$feature_info} containing MZmine IIN IDs
#'   (default: \code{"ion_identities:iin_id"}).
#'   Used when \code{ion_identity = "ion_identity_network"}.
#' @param within_group_dist Distance assigned to pairs within the same ion
#'   identity group (default: \code{0.01}). A small positive value rather than
#'   0 avoids degenerate zero-height clusters in the tree while still pulling
#'   grouped features together. Must be in \code{[0, 1]}.
#' @param save_newick Logical; if \code{TRUE} writes a Newick file
#'   (\code{<outprefix>.nwk}) of the tree (default: \code{FALSE}).
#' @param outprefix File prefix used when \code{save_newick = TRUE}
#'   (default: \code{"feature_tree"}).
#' @return A named list with:
#'   \itemize{
#'     \item \code{hclust}    -- the \code{hclust} object
#'     \item \code{dendrogram} -- the \code{dendrogram} object
#'     \item \code{phylo}     -- the \code{phylo} object (for ape/iTOL)
#'     \item \code{dist_used} -- the (possibly modified) distance matrix used
#'     \item \code{tip_map}   -- data.frame of feature ID -> group assignment
#'                               (\code{NULL} when \code{ion_identity = "none"})
#'   }
#' @export
#' @examplesIf FALSE
#' tree <- FeatureDendrogram(mmo, distance = "dreams")
#' tree_iin <- FeatureDendrogram(mmo, distance = "dreams",
#'                               ion_identity = "ion_identity_network")
# DEPRECATED: superseded by FeatureDendrogram().
# FeatureDendrogram_derep() uses GetDistanceMat() (.dissim slots, passes
# distances directly to hclust). FeatureDendrogram() uses GetSimMat()
# (.sim slots), adds a size guard, and converts S->D (1-S) for hclust.
FeatureDendrogram_derep <- function(
    mmo,
    distance          = "dreams",
    features          = NULL,
    method            = "average",
    ion_identity      = c("none", "correlation", "ion_identity_network"),
    corr_col          = "feature_group",
    iin_col           = "ion_identities:iin_id",
    within_group_dist = 0.01,
    save_newick       = FALSE,
    outprefix         = "feature_tree"
) {
  .require_pkg("ape")
  ion_identity <- match.arg(ion_identity)

  # ------------------------------------------------------------------
  # 1. Retrieve and validate the distance matrix
  # ------------------------------------------------------------------
  mat <- GetDistanceMat(mmo, distance)

  if (!is.matrix(mat))
    stop("Distance matrix '", distance, "' must be a matrix.", call. = FALSE)
  if (is.null(rownames(mat)) || is.null(colnames(mat)))
    stop("Distance matrix must have row and column names (feature IDs).", call. = FALSE)
  if (!identical(rownames(mat), colnames(mat)))
    stop("Distance matrix row names and column names must be identical and in the same order.", call. = FALSE)
  if (any(mat < 0 | mat > 1, na.rm = TRUE))
    stop("Distance matrix contains values outside [0, 1]. ",
         "Ensure dissimilarities are in [0, 1] (e.g. 1 - similarity).", call. = FALSE)

  # ------------------------------------------------------------------
  # 2. Optional feature subsetting
  # ------------------------------------------------------------------
  if (!is.null(features)) {
    keep <- intersect(as.character(features), rownames(mat))
    if (length(keep) < 2)
      stop("Fewer than 2 features remain after subsetting. Check your 'features' vector.", call. = FALSE)
    mat <- mat[keep, keep, drop = FALSE]
  }

  diag(mat) <- 0

  # ------------------------------------------------------------------
  # 3. Ion identity / correlation constraints
  # ------------------------------------------------------------------
  tip_ids <- rownames(mat)
  tip_map <- NULL

  if (ion_identity != "none") {
    if (is.null(mmo$feature_info) || !is.data.frame(mmo$feature_info))
      stop("mmo$feature_info is required when ion_identity != 'none'.", call. = FALSE)
    if (!"id" %in% names(mmo$feature_info))
      stop("mmo$feature_info must contain an 'id' column.", call. = FALSE)

    fi  <- mmo$feature_info
    fi$id <- as.character(fi$id)
    idx <- match(tip_ids, fi$id)

    if (anyNA(idx)) {
      bad <- tip_ids[is.na(idx)]
      stop("Some feature IDs in the distance matrix are missing from mmo$feature_info$id:\n",
           paste(head(bad, 10), collapse = "\n"), call. = FALSE)
    }

    if (ion_identity == "correlation") {
      if (!corr_col %in% names(fi))
        stop("mmo$feature_info is missing column '", corr_col, "' (corr_col).", call. = FALSE)
      grp      <- as.character(fi[[corr_col]][idx])
      grp_type <- "correlation"
    } else {
      if (!iin_col %in% names(fi))
        stop("mmo$feature_info is missing column '", iin_col, "' (iin_col).", call. = FALSE)
      grp      <- as.character(fi[[iin_col]][idx])
      grp_type <- "ion_identity_network"
    }

    grp[is.na(grp) | grp == "" | grp == "NA"] <- NA
    tip_map <- data.frame(id = tip_ids, group_type = grp_type, group_id = grp,
                          stringsAsFactors = FALSE)

    groups <- sort(unique(grp[!is.na(grp)]))
    n_constrained <- 0L
    for (g in groups) {
      ii <- which(grp == g)
      if (length(ii) >= 2) {
        mat[ii, ii] <- within_group_dist
        n_constrained <- n_constrained + length(ii)
      }
    }
    diag(mat) <- 0
    message(sprintf("Ion identity constraint applied: %d groups, %d features affected.",
                    length(groups), n_constrained))
  }

  # ------------------------------------------------------------------
  # 4. Build the tree
  # ------------------------------------------------------------------
  hc   <- stats::hclust(stats::as.dist(mat), method = method)
  dend <- stats::as.dendrogram(hc)
  phy  <- ape::as.phylo(hc)

  if (isTRUE(save_newick)) {
    ape::write.tree(phy, file = paste0(outprefix, ".nwk"))
    message("Newick tree written to ", outprefix, ".nwk")
  }

  return(list(
    hclust     = hc,
    dendrogram = dend,
    phylo      = phy,
    dist_used  = mat,
    tip_map    = tip_map
  ))
}


#' FeatureDendrogram (similarity-based)
#'
#' Identical to \code{FeatureDendrogram_derep()} but retrieves the feature
#' \strong{similarity} matrix via \code{GetSimMat()} (.sim slots) and
#' converts to distance (\code{1 - S}) only inside \code{hclust()}.
#' Adds a size guard: stops with a clear message if n > 10,000 features,
#' because \code{as.dist()} on an n x n matrix is infeasible at that scale.
#' Ion identity constraints are expressed as high similarity
#' (\code{1 - within_group_dist}) rather than low distance.
#'
#' @inheritParams FeatureDendrogram_derep
#' @export
FeatureDendrogram <- function(
    mmo,
    distance          = "dreams",
    features          = NULL,
    method            = "average",
    ion_identity      = c("none", "correlation", "ion_identity_network"),
    corr_col          = "feature_group",
    iin_col           = "ion_identities:iin_id",
    within_group_dist = 0.01,
    save_newick       = FALSE,
    outprefix         = "feature_tree"
) {
  .require_pkg("ape")
  ion_identity <- match.arg(ion_identity)

  # ------------------------------------------------------------------
  # 1. Retrieve similarity matrix and validate
  # ------------------------------------------------------------------
  mat <- GetSimMat(mmo, distance)

  if (!is.matrix(mat) && !inherits(mat, "Matrix"))
    stop("Similarity matrix '", distance, "' must be a matrix.", call. = FALSE)
  if (is.null(rownames(mat)) || is.null(colnames(mat)))
    stop("Similarity matrix must have row and column names (feature IDs).", call. = FALSE)
  if (!identical(rownames(mat), colnames(mat)))
    stop("Similarity matrix row names and column names must be identical and in the same order.", call. = FALSE)

  # Size guard: hclust requires a full dense n x n distance matrix
  n_feat <- nrow(mat)
  if (n_feat > 10000L) {
    stop(
      "FeatureDendrogram requires a dense distance matrix and cannot run at this scale ",
      "(n = ", n_feat, " features).\n",
      "Use filter_mmo() to reduce the feature count below 10,000, then call FeatureDendrogram().",
      call. = FALSE
    )
  }

  # Densify if sparse (safe here — size guard passed)
  mat <- as.matrix(mat)

  if (any(mat < 0 | mat > 1, na.rm = TRUE))
    stop("Similarity matrix contains values outside [0, 1].", call. = FALSE)

  # ------------------------------------------------------------------
  # 2. Optional feature subsetting
  # ------------------------------------------------------------------
  if (!is.null(features)) {
    keep <- intersect(as.character(features), rownames(mat))
    if (length(keep) < 2)
      stop("Fewer than 2 features remain after subsetting. Check your 'features' vector.", call. = FALSE)
    mat <- mat[keep, keep, drop = FALSE]
  }

  diag(mat) <- 1  # self-similarity = 1

  # ------------------------------------------------------------------
  # 3. Ion identity / correlation constraints
  # ------------------------------------------------------------------
  tip_ids <- rownames(mat)
  tip_map <- NULL

  if (ion_identity != "none") {
    if (is.null(mmo$feature_info) || !is.data.frame(mmo$feature_info))
      stop("mmo$feature_info is required when ion_identity != 'none'.", call. = FALSE)
    if (!"id" %in% names(mmo$feature_info))
      stop("mmo$feature_info must contain an 'id' column.", call. = FALSE)

    fi  <- mmo$feature_info
    fi$id <- as.character(fi$id)
    idx <- match(tip_ids, fi$id)

    if (anyNA(idx)) {
      bad <- tip_ids[is.na(idx)]
      stop("Some feature IDs in the similarity matrix are missing from mmo$feature_info$id:\n",
           paste(head(bad, 10), collapse = "\n"), call. = FALSE)
    }

    if (ion_identity == "correlation") {
      if (!corr_col %in% names(fi))
        stop("mmo$feature_info is missing column '", corr_col, "' (corr_col).", call. = FALSE)
      grp      <- as.character(fi[[corr_col]][idx])
      grp_type <- "correlation"
    } else {
      if (!iin_col %in% names(fi))
        stop("mmo$feature_info is missing column '", iin_col, "' (iin_col).", call. = FALSE)
      grp      <- as.character(fi[[iin_col]][idx])
      grp_type <- "ion_identity_network"
    }

    grp[is.na(grp) | grp == "" | grp == "NA"] <- NA
    tip_map <- data.frame(id = tip_ids, group_type = grp_type, group_id = grp,
                          stringsAsFactors = FALSE)

    groups <- sort(unique(grp[!is.na(grp)]))
    n_constrained <- 0L
    for (g in groups) {
      ii <- which(grp == g)
      if (length(ii) >= 2) {
        # Force high similarity (= low distance) for ions in the same group
        mat[ii, ii] <- 1 - within_group_dist
        n_constrained <- n_constrained + length(ii)
      }
    }
    diag(mat) <- 1  # restore self-similarity after constraint
    message(sprintf("Ion identity constraint applied: %d groups, %d features affected.",
                    length(groups), n_constrained))
  }

  # ------------------------------------------------------------------
  # 4. Build the tree (convert similarity to distance for hclust)
  # ------------------------------------------------------------------
  dist_mat <- 1 - mat  # S -> D; densification already done above
  hc   <- stats::hclust(stats::as.dist(dist_mat), method = method)
  dend <- stats::as.dendrogram(hc)
  phy  <- ape::as.phylo(hc)

  if (isTRUE(save_newick)) {
    ape::write.tree(phy, file = paste0(outprefix, ".nwk"))
    message("Newick tree written to ", outprefix, ".nwk")
  }

  return(list(
    hclust     = hc,
    dendrogram = dend,
    phylo      = phy,
    dist_used  = dist_mat,  # distance matrix used for clustering (1 - sim)
    sim_used   = mat,       # original similarity matrix
    tip_map    = tip_map
  ))
}


#' PlotFeatureDendrogram
#'
#' Plot a feature dendrogram produced by \code{FeatureDendrogram()} using
#' \pkg{ggtree}, returning a \code{ggplot}/\code{ggtree} object that can be
#' further customised with additional \code{ggplot2} or \code{ggtree} layers.
#'
#' @details
#' \strong{Design philosophy:}
#' The function returns a \code{ggtree} object so you can keep layering
#' modifications after the call:
#' \preformatted{
#'   p <- PlotFeatureDendrogram(tree, mmo)
#'   p + ggtree::geom_tiplab(size = 1.5) +
#'       ggplot2::theme(legend.position = "bottom")
#' }
#'
#' \strong{Tip point coloring (\code{color_by}):}
#' Because most metabolomics trees have hundreds or thousands of tips,
#' text labels are rarely readable. Instead, tiny colored points are drawn at
#' each tip using \code{ggtree::geom_tippoint()} to indicate compound class.
#' Use \code{show_tip_labels = TRUE} to additionally draw feature ID text (not
#' recommended for large trees).
#'
#' \strong{IIN / correlation group highlighting (\code{highlight_groups}):}
#' When \code{tree$tip_map} is present, each ion identity / correlation group
#' is highlighted using \code{ggtree::geom_hilight()}, which shades the MRCA
#' clade of all members of that group.
#'
#' @param tree Output list from \code{FeatureDendrogram()}.
#' @param mmo \code{mmo}. Used for \code{mmo$sirius_annot} when
#'   \code{color_by} is not \code{NULL}.
#' @param color_by Column name in \code{mmo$sirius_annot} to use for tip point
#'   colors (default: \code{"NPC#pathway"}). Set to \code{NULL} to skip
#'   tip coloring.
#' @param layout Tree layout passed to \code{ggtree::ggtree()}. One of
#'   \code{"circular"} (default), \code{"rectangular"}, \code{"fan"},
#'   \code{"slanted"}, \code{"radial"}, or any other layout supported by
#'   \pkg{ggtree}.
#' @param color_branches Logical; whether to color the tree branches
#'   themselves by the annotation in \code{color_by}, in addition to the tip
#'   points (default: \code{FALSE}). Branch color is determined by a
#'   majority-vote of the tip annotations descending from each node, so
#'   mixed-class clades are colored by whichever class is most abundant.
#'   Unannotated or mixed branches fall back to \code{na_color}.
#' @param highlight_groups Logical; whether to highlight ion identity /
#'   correlation groups using \code{ggtree::geom_hilight()} (default:
#'   \code{TRUE}). Only has an effect when \code{tree$tip_map} is not
#'   \code{NULL}.
#' @param palette Qualitative palette name for
#'   \code{colorspace::qualitative_hcl} (default: \code{"Dark 3"}).
#' @param na_color Color for tips/branches with no annotation (default:
#'   \code{"grey70"}).
#' @param tip_size Point size for tip points (default: \code{1}).
#' @param branch_size Line width for tree branches (default: \code{0.5}).
#'   Increase for small trees; decrease for large ones.
#' @param show_tip_labels Logical; draw feature ID text at tips (default:
#'   \code{FALSE}). Only useful for small trees.
#' @param open_angle For \code{layout = "fan"}, the opening angle in degrees
#'   (default: \code{10}).
#' @return A \code{ggtree}/\code{ggplot} object. Add further layers with
#'   \code{+}.
#' @export
#' @examplesIf FALSE
#' tree <- FeatureDendrogram(mmo, distance = "dreams")
#'
#' # Basic circular tree coloured by NPC pathway (tip points only)
#' PlotFeatureDendrogram(tree, mmo)
#'
#' # Color branches too, thicker lines for a small tree
#' PlotFeatureDendrogram(tree, mmo, color_branches = TRUE, branch_size = 1.2)
#'
#' # Rectangular layout, no coloring
#' PlotFeatureDendrogram(tree, mmo, layout = "rectangular", color_by = NULL)
#'
#' # Returned object is a ggplot -- keep layering
#' p <- PlotFeatureDendrogram(tree, mmo)
#' p + ggtree::geom_tiplab(size = 1, align = TRUE)
PlotFeatureDendrogram <- function(
    tree,
    mmo,
    color_by         = "NPC#pathway",
    layout           = "circular",
    color_branches   = FALSE,
    highlight_groups = TRUE,
    palette          = "Dark 3",
    na_color         = "grey70",
    tip_size         = 1,
    branch_size      = 0.5,
    show_tip_labels  = FALSE,
    open_angle       = 10
) {
  .require_pkg("ggtree")
  .require_pkg("ggplot2")
  .require_pkg("colorspace")
  .require_pkg("ape")

  if (!all(c("hclust", "dendrogram", "phylo", "dist_used") %in% names(tree)))
    stop("'tree' must be the list returned by FeatureDendrogram().", call. = FALSE)

  phy     <- tree$phylo
  tip_ids <- phy$tip.label
  n_tips  <- length(tip_ids)

  # ------------------------------------------------------------------
  # 1. Build tip annotation data.frame joined to phylo
  # ------------------------------------------------------------------
  tip_df <- data.frame(
    label       = tip_ids,
    annot_value = NA_character_,
    stringsAsFactors = FALSE
  )

  levels_present <- character(0)

  if (!is.null(color_by)) {
    if (is.null(mmo$sirius_annot) || !is.data.frame(mmo$sirius_annot))
      stop("mmo$sirius_annot not found. Run AddSiriusAnnot() first, or set color_by = NULL.", call. = FALSE)
    if (!"id" %in% names(mmo$sirius_annot))
      stop("mmo$sirius_annot must contain an 'id' column.", call. = FALSE)
    if (!color_by %in% names(mmo$sirius_annot))
      stop("Column '", color_by, "' not found in mmo$sirius_annot.", call. = FALSE)

    annot    <- mmo$sirius_annot
    annot$id <- as.character(annot$id)
    idx      <- match(tip_ids, annot$id)
    raw_vals <- ifelse(is.na(idx), NA_character_,
                       as.character(annot[[color_by]][idx]))
    raw_vals[is.na(raw_vals) | raw_vals == "" | raw_vals == "NA"] <- NA_character_

    tip_df$annot_value  <- raw_vals
    levels_present      <- sort(unique(raw_vals[!is.na(raw_vals)]))
  }

  # ------------------------------------------------------------------
  # 2. Build color scale
  # ------------------------------------------------------------------
  n_levels <- length(levels_present)

  if (n_levels > 0) {
    hcl_cols <- colorspace::qualitative_hcl(n_levels, palette = palette)
    col_vals <- c(stats::setNames(hcl_cols, levels_present), "NA" = na_color)
    # replace NA with a sentinel string for ggplot discrete scale
    tip_df$annot_value[is.na(tip_df$annot_value)] <- "No annotation"
    col_vals["No annotation"] <- na_color
    col_vals <- col_vals[names(col_vals) != "NA"]
  }

  # ------------------------------------------------------------------
  # 3. Build node-level annotation for branch coloring (majority vote)
  # ------------------------------------------------------------------
  node_annot_df <- NULL

  if (!is.null(color_by) && n_levels > 0 && isTRUE(color_branches)) {
    n_nodes  <- phy$Nnode
    n_total  <- n_tips + n_nodes
    node_class <- rep(NA_character_, n_total)

    # seed with tip values (tips are 1:n_tips in ape ordering)
    tip_order <- match(phy$tip.label, tip_df$label)
    node_class[seq_len(n_tips)] <- tip_df$annot_value[tip_order]

    # majority-vote propagation bottom-up through edge matrix
    # ape edge matrix: col1 = parent, col2 = child
    edge_mat <- phy$edge
    # process in reverse edge order (leaves first)
    for (i in rev(seq_len(nrow(edge_mat)))) {
      parent <- edge_mat[i, 1]
      child  <- edge_mat[i, 2]
      if (is.na(node_class[parent])) {
        # look at all children of this parent
        children <- edge_mat[edge_mat[, 1] == parent, 2]
        child_vals <- node_class[children]
        child_vals <- child_vals[!is.na(child_vals) & child_vals != "No annotation"]
        if (length(child_vals) > 0) {
          tbl <- table(child_vals)
          node_class[parent] <- names(tbl)[which.max(tbl)]
        }
      }
    }

    # replace remaining NA with sentinel
    node_class[is.na(node_class)] <- "No annotation"

    node_annot_df <- data.frame(
      node        = seq_len(n_total),
      branch_class = node_class,
      stringsAsFactors = FALSE
    )
  }

  # ------------------------------------------------------------------
  # 4. Base ggtree
  # ------------------------------------------------------------------
  if (layout == "fan") {
    p <- ggtree::ggtree(phy, layout = "fan", open.angle = open_angle,
                        size = branch_size)
  } else {
    p <- ggtree::ggtree(phy, layout = layout, size = branch_size)
  }

  # ------------------------------------------------------------------
  # 5. Branch coloring (attach node annotation and map to color aes)
  # ------------------------------------------------------------------
  if (!is.null(color_by) && n_levels > 0 && isTRUE(color_branches) &&
      !is.null(node_annot_df)) {
    p <- ggtree::`%<+%`(p, node_annot_df) +
      ggplot2::aes(color = .data$branch_class) +
      ggplot2::scale_color_manual(
        name     = color_by,
        values   = col_vals,
        na.value = na_color
      )
  }

  # ------------------------------------------------------------------
  # 6. Attach tip annotation and add tip points
  # ------------------------------------------------------------------
  if (!is.null(color_by) && n_levels > 0) {
    if (isTRUE(color_branches)) {
      # color scale already set; just add tip points using same aesthetic
      p <- ggtree::`%<+%`(p, tip_df) +
        ggtree::geom_tippoint(
          ggplot2::aes(color = .data$annot_value),
          size  = tip_size,
          na.rm = TRUE
        )
    } else {
      # tip points only -- set color scale here
      p <- ggtree::`%<+%`(p, tip_df) +
        ggtree::geom_tippoint(
          ggplot2::aes(color = .data$annot_value),
          size  = tip_size,
          na.rm = TRUE
        ) +
        ggplot2::scale_color_manual(
          name     = color_by,
          values   = col_vals,
          na.value = na_color
        )
    }
  } else {
    p <- p +
      ggtree::geom_tippoint(size = tip_size, color = na_color)
  }

  # ------------------------------------------------------------------
  # 7. Tip labels (optional)
  # ------------------------------------------------------------------
  if (isTRUE(show_tip_labels)) {
    p <- p + ggtree::geom_tiplab(size = 1.5, align = FALSE)
  }

  # ------------------------------------------------------------------
  # 8. IIN / correlation group highlighting via geom_hilight
  # ------------------------------------------------------------------
  if (isTRUE(highlight_groups) && !is.null(tree$tip_map)) {
    tm      <- tree$tip_map
    tm$id   <- as.character(tm$id)
    grp_ids <- sort(unique(tm$group_id[!is.na(tm$group_id)]))

    if (length(grp_ids) > 0) {
      grp_pal <- stats::setNames(
        colorspace::qualitative_hcl(length(grp_ids), palette = "Set 2"),
        grp_ids
      )

      for (g in grp_ids) {
        members <- tm$id[!is.na(tm$group_id) & tm$group_id == g]
        if (length(members) < 2) next
        # find MRCA node in the phylo object
        node_ids <- which(phy$tip.label %in% members)
        if (length(node_ids) < 2) next
        mrca_node <- ape::getMRCA(phy, node_ids)
        if (is.null(mrca_node)) next
        p <- p + ggtree::geom_hilight(
          node  = mrca_node,
          fill  = grp_pal[[g]],
          alpha = 0.25
        )
      }
    }
  }

  # ------------------------------------------------------------------
  # 9. Theme
  # ------------------------------------------------------------------
  p <- p + ggplot2::theme(
    legend.position = "right",
    legend.key.size = ggplot2::unit(0.4, "cm"),
    legend.text     = ggplot2::element_text(size = 7)
  )

  p
}


# NOTE: ExportITOL has not yet been updated to the similarity-based storage
# approach (.sim slots / GetSimMat). It currently uses GetDistanceMat() and
# .dissim slots. Call AddChemDist() before using this function until it is
# updated in a future release.
#' ExportITOL
#'
#' Export a feature dendrogram and companion annotation files for visualisation
#' in iTOL (Interactive Tree of Life). Produces three files:
#' \enumerate{
#'   \item A Newick tree (\code{<outprefix>.nwk})
#'   \item A colour-strip annotation coloured by NPC pathway (or any
#'         \code{mmo$sirius_annot} column) (\code{<outprefix>_colorstrip.txt})
#'   \item A bar-chart annotation showing per-feature prevalence across samples
#'         (\code{<outprefix>_barplot.txt})
#' }
#'
#' @details
#' \strong{How to load in iTOL:}
#' Upload the Newick file at \url{https://itol.embl.de/upload.cgi}, then drag
#' and drop the two annotation files onto the displayed tree. Select
#' "Circular" layout and enable "Display" for each dataset in the
#' Datasets panel.
#'
#' \strong{Prevalence calculation:}
#' Feature prevalence is the proportion of samples in which a feature is
#' detected (abundance > 0). Features absent from the PA table are assigned
#' prevalence 0. This value drives the bar-chart length in iTOL.
#'
#' \strong{Colour-strip annotation:}
#' Each tip is assigned the colour of its NPC pathway (or other chosen column)
#' from \code{mmo$sirius_annot}. Features with no annotation are coloured
#' \code{na_color}. The legend is written into the file header so iTOL renders
#' it automatically.
#'
#' @param tree Output list from \code{FeatureDendrogram()}.
#' @param mmo \code{mmo}. Must contain \code{mmo$sirius_annot} (for colour
#'   strip) and \code{mmo$feature_presence} or \code{mmo$feature_data} (for
#'   prevalence bar chart).
#' @param outprefix File path prefix for output files (default:
#'   \code{"itol_export"}). Directories in the path must exist.
#' @param color_by Column name in \code{mmo$sirius_annot} used for the colour
#'   strip (default: \code{"NPC#pathway"}).
#' @param palette Qualitative palette name for
#'   \code{colorspace::qualitative_hcl} (default: \code{"Dark 3"}).
#' @param na_color Hex colour for features with no annotation
#'   (default: \code{"#CCCCCC"}).
#' @return Invisibly returns a named character vector of written file paths.
#' @export
#' @examplesIf FALSE
#' tree <- FeatureDendrogram(mmo, distance = "dreams")
#' ExportITOL(tree, mmo, outprefix = "output/itol_dreams")
#' # Upload output/itol_dreams.nwk to iTOL, then drag and drop the two .txt files
ExportITOL <- function(
    tree,
    mmo,
    outprefix = "itol_export",
    color_by  = "NPC#pathway",
    palette   = "Dark 3",
    na_color  = "#CCCCCC"
) {
  .require_pkg("ape")
  .require_pkg("colorspace")

  if (!all(c("phylo", "dist_used") %in% names(tree)))
    stop("'tree' must be the list returned by FeatureDendrogram().", call. = FALSE)

  phy     <- tree$phylo
  tip_ids <- phy$tip.label

  # ------------------------------------------------------------------
  # 1. Newick file
  # ------------------------------------------------------------------
  nwk_path <- paste0(outprefix, ".nwk")
  ape::write.tree(phy, file = nwk_path)
  message("Tree written to: ", nwk_path)

  # ------------------------------------------------------------------
  # 2. Colour-strip annotation (DATASET_COLORSTRIP)
  # ------------------------------------------------------------------
  if (is.null(mmo$sirius_annot) || !is.data.frame(mmo$sirius_annot))
    stop("mmo$sirius_annot not found. Run AddSiriusAnnot() first.", call. = FALSE)
  if (!color_by %in% names(mmo$sirius_annot))
    stop("Column '", color_by, "' not found in mmo$sirius_annot.", call. = FALSE)

  annot      <- mmo$sirius_annot
  annot$id   <- as.character(annot$id)
  idx        <- match(tip_ids, annot$id)
  raw_vals   <- ifelse(is.na(idx), NA_character_, as.character(annot[[color_by]][idx]))
  raw_vals[is.na(raw_vals) | raw_vals == "" | raw_vals == "NA"] <- NA_character_

  levels_present <- sort(unique(raw_vals[!is.na(raw_vals)]))
  if (length(levels_present) == 0)
    warning("No features had a non-missing value for '", color_by,
            "'. All tips will be coloured na_color.", call. = FALSE)

  pal <- if (length(levels_present) > 0) {
    stats::setNames(
      colorspace::qualitative_hcl(length(levels_present), palette = palette),
      levels_present
    )
  } else {
    character(0)
  }
  header_color  <- if (length(pal) > 0) pal[[1L]] else na_color
  tip_colors    <- ifelse(is.na(raw_vals), na_color, unname(pal[raw_vals]))
  tip_labels    <- ifelse(is.na(raw_vals), "No annotation", raw_vals)

  legend_colors <- paste(c(unname(pal), na_color), collapse = "\t")
  legend_labels <- paste(c(levels_present, "No annotation"), collapse = "\t")
  legend_shapes <- paste(rep("1", length(levels_present) + 1L), collapse = "\t")

  strip_path <- paste0(outprefix, "_colorstrip.txt")
  strip_lines <- c(
    "DATASET_COLORSTRIP",
    "SEPARATOR TAB",
    paste0("DATASET_LABEL\t", color_by),
    paste0("COLOR\t", header_color),
    "STRIP_WIDTH\t30",
    "SHOW_STRIP_LABELS\t0",
    paste0("LEGEND_TITLE\t", color_by),
    paste0("LEGEND_SHAPES\t", legend_shapes),
    paste0("LEGEND_COLORS\t", legend_colors),
    paste0("LEGEND_LABELS\t", legend_labels),
    "DATA",
    paste(tip_ids, tip_colors, tip_labels, sep = "\t")
  )
  writeLines(strip_lines, strip_path)
  message("Colour-strip annotation written to: ", strip_path)

  # ------------------------------------------------------------------
  # 3. Prevalence bar chart (DATASET_SIMPLEBAR)
  # ------------------------------------------------------------------
  if (!is.null(mmo$feature_presence) && is.data.frame(mmo$feature_presence)) {
    pa <- mmo$feature_presence
  } else if (!is.null(mmo$feature_data) && is.data.frame(mmo$feature_data)) {
    pa <- mmo$feature_data
  } else {
    stop("mmo must contain feature_presence or feature_data to compute prevalence.", call. = FALSE)
  }

  pa_ids   <- as.character(pa$id)
  samp_mat <- as.matrix(pa[, -c(1, 2), drop = FALSE])
  n_samp   <- ncol(samp_mat)
  detected <- rowSums(samp_mat > 0, na.rm = TRUE)
  prev_map <- stats::setNames(detected / n_samp, pa_ids)

  prevalence <- ifelse(tip_ids %in% names(prev_map), prev_map[tip_ids], 0)
  prevalence <- round(as.numeric(prevalence), 4)

  bar_path <- paste0(outprefix, "_barplot.txt")
  bar_lines <- c(
    "DATASET_SIMPLEBAR",
    "SEPARATOR TAB",
    "DATASET_LABEL\tPrevalence",
    "COLOR\t#444444",
    "FIELD_LABELS\tPrevalence (prop. samples detected)",
    "DATA",
    paste(tip_ids, prevalence, sep = "\t")
  )
  writeLines(bar_lines, bar_path)
  message("Bar-chart annotation written to: ", bar_path)

  invisible(c(tree = nwk_path, colorstrip = strip_path, barplot = bar_path))
}


#' ExportCytoscape
#'
#' Export a \code{mmo} as a node table and edge table for visualisation in
#' Cytoscape. Edges are derived from a stored feature dissimilarity matrix.
#' Nodes carry abundance statistics and any annotation columns present in
#' \code{mmo$sirius_annot} and \code{mmo$feature_info}.
#'
#' @details
#' \strong{Edge filtering:}
#' Two complementary filters control which edges are retained. Both are
#' applied together (a pair must pass both to be included):
#' \itemize{
#'   \item \code{sim_threshold} -- drop all edges with similarity below this
#'     value. Higher values produce sparser, higher-confidence networks.
#'   \item \code{top_k} -- for each node, retain only its \code{top_k} most
#'     similar neighbours (by similarity). \code{NULL} keeps all edges that
#'     pass the threshold. This is useful for preventing highly-connected hub
#'     nodes from dominating the layout.
#' }
#' A warning is issued when the retained edge count exceeds 50 000, as large
#' networks can be slow to render in Cytoscape.
#'
#' \strong{Node table columns:}
#' The node table always includes \code{id}, \code{prevalence} (proportion of
#' samples with detected abundance), and one \code{mean_<group>} column per
#' group. Any columns present in \code{mmo$feature_info} and
#' \code{mmo$sirius_annot} are appended automatically -- the function does not
#' assume specific column names. Column names are sanitised to replace
#' characters that cause problems in Cytoscape (spaces, \code{#}, \code{:}).
#'
#' \strong{Loading in Cytoscape:}
#' \enumerate{
#'   \item File -> Import -> Network from File -> select \code{_edges.csv}.
#'     Map \code{source} as Source Node, \code{target} as Target Node,
#'     \code{similarity} as Edge Attribute.
#'   \item File -> Import -> Table from File -> select \code{_nodes.csv}.
#'     Map \code{id} as the Key column matching node names.
#'   \item In the Style panel, map \code{Fill Color} to any annotation column
#'     (e.g. \code{NPC_pathway}) to color nodes by compound class.
#' }
#'
#' @param mmo \code{mmo}. Must contain \code{mmo$feature_data} and
#'   \code{mmo$metadata}. \code{mmo$feature_info} and
#'   \code{mmo$sirius_annot} are used when present.
#' @param distance Name of the dissimilarity matrix to use for edges
#'   (default: \code{"dreams"}). Passed to \code{GetDistanceMat()}.
#' @param outprefix File path prefix for output CSVs
#'   (default: \code{"cytoscape_export"}).
#' @param sim_threshold Minimum pairwise similarity (= 1 - dissimilarity) for
#'   an edge to be retained (default: \code{0.7}). Set to \code{0} to disable.
#' @param top_k Integer or \code{NULL}. If not \code{NULL}, each node retains
#'   only its \code{top_k} most similar neighbours after threshold filtering
#'   (default: \code{NULL} = no k limit).
#' @param group_col Metadata column defining groups for per-group mean
#'   abundance columns (default: \code{"group"}).
#' @param sample_col Metadata column matching sample IDs to feature columns
#'   (default: \code{"sample"}).
#' @return Invisibly returns a named list:
#'   \code{nodes} (data.frame), \code{edges} (data.frame),
#'   \code{nodes_path}, \code{edges_path} (file paths).
#' @export
#' @examplesIf FALSE
#' # Default: similarity >= 0.7
#' ExportCytoscape(mmo, distance = "dreams",
#'                 outprefix = "output/network")
#'
#' # Sparser network: top 5 neighbours per node, similarity >= 0.6
#' ExportCytoscape(mmo, distance = "dreams",
#'                 outprefix = "output/network_k5",
#'                 sim_threshold = 0.6, top_k = 5)
# NOTE: ExportCytoscape has not yet been updated to the similarity-based
# storage approach (.sim slots / GetSimMat). It currently uses GetDistanceMat()
# and .dissim slots. Call AddChemDist() before using this function until it is
# updated in a future release.
ExportCytoscape <- function(
    mmo,
    distance      = "dreams",
    outprefix     = "cytoscape_export",
    sim_threshold = 0.7,
    top_k         = NULL,
    group_col     = "group",
    sample_col    = "sample"
) {
  if (is.null(mmo$feature_data) || !is.data.frame(mmo$feature_data))
    stop("mmo$feature_data is required.", call. = FALSE)
  if (is.null(mmo$metadata) || !is.data.frame(mmo$metadata))
    stop("mmo$metadata is required.", call. = FALSE)

  feat_ids <- as.character(mmo$feature_data$id)

  # ------------------------------------------------------------------
  # 1. Edge table
  # ------------------------------------------------------------------
  dist_mat  <- GetDistanceMat(mmo, distance)
  keep      <- intersect(feat_ids, rownames(dist_mat))
  dist_mat  <- dist_mat[keep, keep, drop = FALSE]
  n         <- nrow(dist_mat)

  # upper triangle -- all pairs, no self-loops
  ut        <- which(upper.tri(dist_mat), arr.ind = TRUE)
  sim_vals  <- 1 - dist_mat[ut]

  # filter 1: similarity threshold
  pass      <- sim_vals >= sim_threshold
  src       <- rownames(dist_mat)[ut[, 1]]
  tgt       <- rownames(dist_mat)[ut[, 2]]
  sim_pass  <- sim_vals[pass]
  src_pass  <- src[pass]
  tgt_pass  <- tgt[pass]

  # filter 2: top-k per node (keep top_k neighbours per node by similarity)
  if (!is.null(top_k)) {
    top_k <- as.integer(top_k)
    if (top_k < 1L) stop("top_k must be >= 1.", call. = FALSE)
    # for each node, rank its edges by similarity descending
    # an edge is kept if it is in the top_k for either endpoint
    rank_src <- stats::ave(-sim_pass, src_pass, FUN = function(x) rank(x, ties.method = "first"))
    rank_tgt <- stats::ave(-sim_pass, tgt_pass, FUN = function(x) rank(x, ties.method = "first"))
    in_topk  <- (rank_src <= top_k) | (rank_tgt <= top_k)
    sim_pass <- sim_pass[in_topk]
    src_pass <- src_pass[in_topk]
    tgt_pass <- tgt_pass[in_topk]
  }

  n_edges <- length(sim_pass)
  if (n_edges > 50000L)
    warning(sprintf(
      "Edge table contains %d edges. Networks this large can be slow in Cytoscape. ",
      "Consider raising sim_threshold or lowering top_k.", n_edges
    ), call. = FALSE)

  message(sprintf(
    "Edge table: %d edges retained (sim >= %.2f%s, from %d feature pairs).",
    n_edges, sim_threshold,
    if (!is.null(top_k)) paste0(", top_k = ", top_k) else "",
    n * (n - 1L) / 2L
  ))

  edges <- data.frame(
    source          = src_pass,
    target          = tgt_pass,
    similarity      = round(sim_pass, 4),
    distance_method = distance,
    stringsAsFactors = FALSE
  )

  # ------------------------------------------------------------------
  # 2. Node table
  # ------------------------------------------------------------------
  .sanitise_colname <- function(x) {
    x <- gsub("[#: ]+", "_", x)
    x <- gsub("_+$", "", x)
    x
  }

  nodes <- data.frame(id = feat_ids, stringsAsFactors = FALSE)

  # feature_info -- append all columns except id, sanitise names
  if (!is.null(mmo$feature_info) && is.data.frame(mmo$feature_info)) {
    fi      <- mmo$feature_info
    fi$id   <- as.character(fi$id)
    fi_sub  <- fi[match(feat_ids, fi$id), , drop = FALSE]
    fi_sub  <- fi_sub[, setdiff(names(fi_sub), "id"), drop = FALSE]
    names(fi_sub) <- .sanitise_colname(names(fi_sub))
    # drop columns already in nodes to avoid duplicates
    fi_sub  <- fi_sub[, setdiff(names(fi_sub), names(nodes)), drop = FALSE]
    nodes   <- cbind(nodes, fi_sub)
  }

  # prevalence from feature_data
  samp_mat <- as.matrix(mmo$feature_data[, -c(1, 2), drop = FALSE])
  rownames(samp_mat) <- feat_ids
  n_samp   <- ncol(samp_mat)
  nodes$prevalence <- round(rowSums(samp_mat > 0, na.rm = TRUE) / n_samp, 4)

  # per-group mean abundance
  md     <- mmo$metadata
  groups <- sort(unique(as.character(md[[group_col]])))
  for (g in groups) {
    samps_g  <- as.character(md[[sample_col]][md[[group_col]] == g])
    samps_g  <- intersect(samps_g, colnames(samp_mat))
    col_name <- paste0("mean_", .sanitise_colname(g))
    nodes[[col_name]] <- if (length(samps_g) > 0)
      round(rowMeans(samp_mat[, samps_g, drop = FALSE], na.rm = TRUE), 2)
    else
      NA_real_
  }

  # sirius_annot -- append all columns except id, sanitise names
  if (!is.null(mmo$sirius_annot) && is.data.frame(mmo$sirius_annot)) {
    annot      <- mmo$sirius_annot
    annot$id   <- as.character(annot$id)
    idx        <- match(feat_ids, annot$id)
    annot_sub  <- annot[idx, , drop = FALSE]
    annot_sub  <- annot_sub[, setdiff(names(annot_sub), "id"), drop = FALSE]
    names(annot_sub) <- .sanitise_colname(names(annot_sub))
    annot_sub  <- annot_sub[, setdiff(names(annot_sub), names(nodes)), drop = FALSE]
    # blank -> NA
    for (col in names(annot_sub)) {
      v <- as.character(annot_sub[[col]])
      v[v == "" | v == "NA"] <- NA_character_
      annot_sub[[col]] <- v
    }
    nodes <- cbind(nodes, annot_sub)
  }

  # ------------------------------------------------------------------
  # 3. Write files
  # ------------------------------------------------------------------
  nodes_path <- paste0(outprefix, "_nodes.csv")
  edges_path <- paste0(outprefix, "_edges.csv")

  utils::write.csv(nodes, nodes_path, row.names = FALSE, na = "")
  message("Node table: ", nrow(nodes), " nodes, ", ncol(nodes),
          " columns -> ", nodes_path)

  utils::write.csv(edges, edges_path, row.names = FALSE, na = "")
  message("Edge table: ", nrow(edges), " edges -> ", edges_path)

  invisible(list(nodes      = nodes,      edges      = edges,
                 nodes_path = nodes_path, edges_path = edges_path))
}
