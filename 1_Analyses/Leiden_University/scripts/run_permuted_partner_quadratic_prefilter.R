#!/usr/bin/env Rscript

# ============================================================================
# PARC empirical-null prefilter analysis
# Partners: Aristotle_University, Ghent_University, BPI
#
# Purpose:
#   For each of 100 permuted metadata files:
#     1. Match metadata to the unchanged raw count matrix.
#     2. Apply each partner's original preprocessing/normalization choice.
#     3. Run DRomics itemselect(..., select.method = "quadratic", FDR = 0.01).
#     4. Save the complete quadratic-test object and compact result tables.
#
# Important:
#   - No samples are removed.
#   - Concentration labels come from each permuted metadata file.
#   - Expression values/counts are never permuted.
#   - Partner-specific choices are kept separate below.
#   - Adapt only the paths and column-name candidates in USER SETTINGS.
# ============================================================================


# ---- Packages ----------------------------------------------------------------

required_packages <- c("DRomics", "edgeR", "DESeq2")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(DRomics)
  library(edgeR)
  library(DESeq2)
})


# ---- USER SETTINGS ------------------------------------------------------------

# Raw count file shared by all analyses.
# Rows should be probes/genes and columns should be samples.
COUNT_FILE <- "data/EUT046/EUT046_counts_per_gene_per_sample_raw_final.csv"

# Folder containing the 100 separately permuted metadata files.
PERMUTED_METADATA_DIR <- "output/EUT046/permuted_metadata"

# Output root.
OUTPUT_DIR <- "output/EUT046/permuted_quadratic_prefilter"

# File patterns.
METADATA_PATTERN <- "\\.(csv|tsv|txt)$"

# Separator is detected from the extension:
#   .csv      -> comma
#   .tsv/.txt -> tab

# Candidate metadata columns. The first matching name is used.
SAMPLE_COLUMN_CANDIDATES <- "SAMPLE_ID"

DOSE_COLUMN_CANDIDATES <- "CONCENTRATION"

TIME_COLUMN_CANDIDATES <- "TIME"

# Raw-count identifier column candidates.
COUNT_ID_COLUMN_CANDIDATES <- "Column1"

# Optional map from probe IDs to the identifier used by Aristotle.
# Aristotle's shared code joins probes to Entrez IDs and sums duplicated Entrez IDs.
#
# Set to NULL only when row names in COUNT_FILE already contain the final IDs.
ARISTOTLE_ANNOTATION_FILE <- "data/Manifests/Temposeq_manifest_Human_WT_1.2_realignment_96percent_2023-05-15.txt"
ARISTOTLE_PROBE_COLUMN_CANDIDATES <- c("Probe_ID", "probe_id", "ProbeID")
ARISTOTLE_GENE_COLUMN_CANDIDATES <- c("entrez_id", "ENTREZID", "EntrezID", "gene_id")

# Transfo method Ghent
GHENT_TRANSFO_METHOD <- "vst"

# Quadratic prefilter threshold used in the partner code.
PREFILTER_FDR <- 0.01

# Set TRUE to overwrite completed runs.
OVERWRITE <- FALSE


# ---- General helpers ----------------------------------------------------------

first_existing_name <- function(nms, candidates, label) {
  hit <- candidates[candidates %in% nms]
  if (length(hit) == 0L) {
    stop(
      "Could not identify ", label, ". Available columns: ",
      paste(nms, collapse = ", ")
    )
  }
  hit[[1L]]
}


read_table_auto <- function(path, check.names = FALSE) {
  ext <- tolower(tools::file_ext(path))
  
  if (ext == "csv") {
    read.csv(
      path,
      check.names = check.names,
      stringsAsFactors = FALSE
    )
  } else {
    read.delim(
      path,
      check.names = check.names,
      stringsAsFactors = FALSE
    )
  }
}


clean_numeric_dose <- function(x) {
  if (is.numeric(x)) {
    return(as.numeric(x))
  }
  
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", "null", "NULL")] <- NA_character_
  
  # Preserve decimal points and scientific notation while removing units.
  cleaned <- gsub(",", ".", x, fixed = TRUE)
  cleaned <- gsub("[^0-9eE+.\\-]", "", cleaned, perl = TRUE)
  
  out <- suppressWarnings(as.numeric(cleaned))
  
  if (anyNA(out)) {
    bad <- unique(x[is.na(out)])
    stop(
      "Dose conversion failed for: ",
      paste(bad, collapse = ", ")
    )
  }
  
  out
}


standardize_time <- function(x) {
  trimws(as.character(x))
}


read_count_matrix <- function(path) {
  if (!file.exists(path)) {
    stop("COUNT_FILE does not exist: ", path)
  }
  
  x <- read_table_auto(path, check.names = FALSE)
  
  id_col <- first_existing_name(
    names(x),
    COUNT_ID_COLUMN_CANDIDATES,
    "the count identifier column"
  )
  
  ids <- as.character(x[[id_col]])
  x[[id_col]] <- NULL
  
  # Convert every sample column to numeric without changing its order.
  x[] <- lapply(x, function(z) {
    out <- suppressWarnings(as.numeric(as.character(z)))
    if (anyNA(out) && !all(is.na(z))) {
      stop("Non-numeric values found in count column.")
    }
    out
  })
  
  mat <- as.matrix(x)
  rownames(mat) <- ids
  
  if (anyDuplicated(colnames(mat))) {
    stop("The count matrix contains duplicated sample column names.")
  }
  
  mat
}


read_permuted_metadata <- function(path) {
  meta <- read_table_auto(path, check.names = FALSE)
  
  sample_col <- first_existing_name(
    names(meta),
    SAMPLE_COLUMN_CANDIDATES,
    "the sample column in metadata"
  )
  dose_col <- first_existing_name(
    names(meta),
    DOSE_COLUMN_CANDIDATES,
    "the concentration/dose column in metadata"
  )
  time_col <- first_existing_name(
    names(meta),
    TIME_COLUMN_CANDIDATES,
    "the timepoint column in metadata"
  )
  
  out <- data.frame(
    sample = as.character(meta[[sample_col]]),
    dose = clean_numeric_dose(meta[[dose_col]]),
    timepoint = standardize_time(meta[[time_col]]),
    stringsAsFactors = FALSE
  )
  
  if (anyDuplicated(out$sample)) {
    stop("Duplicated sample names in metadata file: ", basename(path))
  }
  
  if (anyNA(out$sample) || any(out$sample == "")) {
    stop("Missing sample names in metadata file: ", basename(path))
  }
  
  if (anyNA(out$dose)) {
    stop("Missing doses in metadata file: ", basename(path))
  }
  
  if (anyNA(out$timepoint) || any(out$timepoint == "")) {
    stop("Missing timepoints in metadata file: ", basename(path))
  }
  
  out
}


align_counts_and_metadata <- function(counts, meta) {
  missing_in_counts <- setdiff(meta$sample, colnames(counts))
  extra_in_counts <- setdiff(colnames(counts), meta$sample)
  
  if (length(missing_in_counts) > 0L) {
    stop(
      "Metadata samples absent from count matrix: ",
      paste(missing_in_counts, collapse = ", ")
    )
  }
  
  # All partner samples are retained. Therefore extra count columns are an error.
  if (length(extra_in_counts) > 0L) {
    stop(
      "Count-matrix samples absent from metadata: ",
      paste(extra_in_counts, collapse = ", ")
    )
  }
  
  counts <- counts[, meta$sample, drop = FALSE]
  
  stopifnot(identical(colnames(counts), meta$sample))
  
  list(counts = counts, metadata = meta)
}


safe_dir_create <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}


selected_item_table <- function(selection) {
  if (length(selection$selectindex) == 0L) {
    return(data.frame(
      item = character(0),
      row_index = integer(0),
      adjusted_p_value = numeric(0)
    ))
  }
  
  item_names <- rownames(selection$omicdata$data)
  if (is.null(item_names)) {
    item_names <- as.character(selection$omicdata$item)
  }
  
  data.frame(
    item = item_names[selection$selectindex],
    row_index = selection$selectindex,
    adjusted_p_value = as.numeric(selection$adjpvalue),
    stringsAsFactors = FALSE
  )
}


all_item_table <- function(selection) {
  # DRomics versions differ slightly in stored object structure.
  # Save the complete object regardless; export all adjusted p-values when present.
  candidate_names <- c(
    "all.adjpvalue", "alladjpvalue", "adjpvalue.all",
    "pvalue.adjusted", "allpvalue"
  )
  
  stored_name <- candidate_names[candidate_names %in% names(selection)][1L]
  
  if (is.na(stored_name)) {
    return(NULL)
  }
  
  values <- selection[[stored_name]]
  data.frame(
    item = names(values),
    adjusted_p_value = as.numeric(values),
    stringsAsFactors = FALSE
  )
}


run_quadratic <- function(omic_object) {
  DRomics::itemselect(
    omic_object,
    select.method = "quadratic",
    FDR = PREFILTER_FDR
  )
}


save_selection <- function(selection, output_prefix, run_info) {
  saveRDS(selection, paste0(output_prefix, "_itemselect.rds"))
  
  write.csv(
    selected_item_table(selection),
    paste0(output_prefix, "_selected_items.csv"),
    row.names = FALSE
  )
  
  all_items <- all_item_table(selection)
  if (!is.null(all_items)) {
    write.csv(
      all_items,
      paste0(output_prefix, "_all_items.csv"),
      row.names = FALSE
    )
  }
  
  run_info$n_selected <- length(selection$selectindex)
  write.csv(
    run_info,
    paste0(output_prefix, "_run_info.csv"),
    row.names = FALSE
  )
}


# ---- Partner-specific preprocessing -------------------------------------------

read_aristotle_annotation <- function() {
  if (is.null(ARISTOTLE_ANNOTATION_FILE)) {
    return(NULL)
  }
  
  if (!file.exists(ARISTOTLE_ANNOTATION_FILE)) {
    stop(
      "Aristotle annotation file does not exist: ",
      ARISTOTLE_ANNOTATION_FILE
    )
  }
  
  ann <- read_table_auto(ARISTOTLE_ANNOTATION_FILE, check.names = FALSE)
  
  probe_col <- first_existing_name(
    names(ann),
    ARISTOTLE_PROBE_COLUMN_CANDIDATES,
    "Aristotle's probe-ID annotation column"
  )
  gene_col <- first_existing_name(
    names(ann),
    ARISTOTLE_GENE_COLUMN_CANDIDATES,
    "Aristotle's final gene-ID annotation column"
  )
  
  data.frame(
    probe_id = as.character(ann[[probe_col]]),
    gene_id = as.character(ann[[gene_col]]),
    stringsAsFactors = FALSE
  )
}


aggregate_by_id <- function(mat, ids) {
  keep <- !is.na(ids) & ids != ""
  mat <- mat[keep, , drop = FALSE]
  ids <- ids[keep]
  
  rowsum(mat, group = ids, reorder = FALSE, na.rm = TRUE)
}


prepare_aristotle <- function(raw_counts, annotation) {
  # Partner code:
  #   calcNormFactors(method = "upperquartile", p = 0.75)
  #   CPM
  #   retain rows with CPM >= 1 in at least ceiling(0.75 * 3) samples
  #   map probes to Entrez ID and sum duplicate Entrez IDs
  #   pseudocount = 10% of smallest non-zero value
  #   log2 transformation
  dge <- edgeR::DGEList(counts = raw_counts)
  dge <- edgeR::calcNormFactors(
    dge,
    method = "upperquartile",
    p = 0.75
  )
  
  cpm_values <- edgeR::cpm(
    dge,
    normalized.lib.sizes = TRUE,
    log = FALSE
  )
  
  min_replicates_threshold <- ceiling(0.75 * 3)
  keep <- rowSums(cpm_values >= 1) >= min_replicates_threshold
  cpm_values <- cpm_values[keep, , drop = FALSE]
  
  if (!is.null(annotation)) {
    gene_id <- annotation$gene_id[
      match(rownames(cpm_values), annotation$probe_id)
    ]
    cpm_values <- aggregate_by_id(cpm_values, gene_id)
  } else if (anyDuplicated(rownames(cpm_values))) {
    cpm_values <- rowsum(
      cpm_values,
      group = rownames(cpm_values),
      reorder = FALSE
    )
  }
  
  non_zero <- cpm_values[cpm_values > 0]
  if (length(non_zero) == 0L) {
    stop("Aristotle preprocessing produced no non-zero values.")
  }
  
  pseudocount <- min(non_zero) * 0.1
  log2(cpm_values + pseudocount)
}


prepare_ghent <- function(raw_counts, metadata) {
  # Ghent uses raw integer counts and requests VST through DRomics.
  if (anyDuplicated(rownames(raw_counts))) {
    raw_counts <- rowsum(
      raw_counts,
      group = rownames(raw_counts),
      reorder = FALSE
    )
  }
  
  storage.mode(raw_counts) <- "integer"
  raw_counts
}


prepare_bpi <- function(raw_counts) {
  # BPI code aggregates duplicated IDs first, then passes raw, unnormalized
  # integer counts to formatdata4DRomics() and RNAseqdata().
  if (anyDuplicated(rownames(raw_counts))) {
    raw_counts <- rowsum(
      raw_counts,
      group = rownames(raw_counts),
      reorder = FALSE
    )
  }
  
  storage.mode(raw_counts) <- "integer"
  raw_counts
}


# ---- Timepoint analysis --------------------------------------------------------

make_dromics_input <- function(signal_matrix, doses) {
  DRomics::formatdata4DRomics(
    signalmatrix = signal_matrix,
    dose = doses,
    samplenames = colnames(signal_matrix)
  )
}


run_partner_timepoint <- function(
    partner,
    signal_matrix,
    doses,
    output_prefix
) {
  formatted <- make_dromics_input(signal_matrix, doses)
  
  if (partner == "Ghent_University") {
    # Explicit VST inside DRomics/DESeq2.
    omic_object <- DRomics::RNAseqdata(
      formatted,
      transfo.method = GHENT_TRANSFO_METHOD
    )
  } else if (partner == "BPI") {
    # BPI omitted transfo.method in the submitted code.
    # DRomics therefore chooses rlog for <30 samples and vst otherwise.
    omic_object <- DRomics::RNAseqdata(formatted)
  } else {
    # Aristotle supplies already transformed continuous data.
    omic_object <- DRomics::continuousomicdata(
      formatted,
      backgrounddose = 0,
      check = TRUE
    )
  }
  
  selection <- run_quadratic(omic_object)
  
  list(
    selection = selection,
    n_selected = length(selection$selectindex)
  )
}


# ---- Main loop ----------------------------------------------------------------

safe_dir_create(OUTPUT_DIR)

metadata_files <- list.files(
  PERMUTED_METADATA_DIR,
  pattern = METADATA_PATTERN,
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(metadata_files) == 0L) {
  stop(
    "No metadata files found in: ",
    PERMUTED_METADATA_DIR
  )
}

metadata_files <- sort(metadata_files)

if (length(metadata_files) != 100L) {
  warning(
    "Expected 100 permuted metadata files, found ",
    length(metadata_files),
    ". All found files will be processed."
  )
}

raw_counts <- read_count_matrix(COUNT_FILE)
aristotle_annotation <- read_aristotle_annotation()

summary_rows <- list()
summary_index <- 1L

for (metadata_file in metadata_files) {
  permutation_id <- tools::file_path_sans_ext(basename(metadata_file))
  message("Permutation: ", permutation_id)
  
  metadata <- read_permuted_metadata(metadata_file)
  aligned <- align_counts_and_metadata(raw_counts, metadata)
  
  counts_aligned <- aligned$counts
  metadata <- aligned$metadata
  
  # Preprocessing is repeated per metadata file so that all objects retain
  # exactly the sample ordering in that permutation's metadata.
  aristotle_matrix <- prepare_aristotle(
    counts_aligned,
    aristotle_annotation
  )
  ghent_counts <- prepare_ghent(
    counts_aligned,
    metadata
  )
  bpi_counts <- prepare_bpi(counts_aligned)
  
  partner_objects <- list(
    Aristotle_University = list(
      matrix = aristotle_matrix,
      method = "upperquartile_CPM_filter_sumIDs_pseudocount_log2"
    ),
    Ghent_University = list(
      matrix = ghent_counts,
      method = paste0("raw_counts_then_DRomics_", GHENT_TRANSFO_METHOD)
    ),
    BPI = list(
      matrix = bpi_counts,
      method = "raw_counts_then_DRomics_RNAseqdata"
    )
  )
  
  for (partner in names(partner_objects)) {
    partner_dir <- file.path(OUTPUT_DIR, partner, permutation_id)
    safe_dir_create(partner_dir)
    
    partner_matrix <- partner_objects[[partner]]$matrix
    
    for (timepoint in unique(metadata$timepoint)) {
      idx <- which(metadata$timepoint == timepoint)
      sample_names <- metadata$sample[idx]
      doses <- metadata$dose[idx]
      
      signal_matrix <- partner_matrix[, sample_names, drop = FALSE]
      
      # DRomics quadratic testing requires multiple dose ranks.
      if (length(unique(doses)) < 2L) {
        warning(
          partner, " / ", permutation_id, " / ", timepoint,
          ": fewer than two unique doses; skipped."
        )
        next
      }
      
      safe_timepoint <- gsub("[^A-Za-z0-9._-]", "_", timepoint)
      output_prefix <- file.path(
        partner_dir,
        paste0("quadratic_", safe_timepoint)
      )
      
      completion_file <- paste0(output_prefix, "_itemselect.rds")
      if (file.exists(completion_file) && !OVERWRITE) {
        message("  Existing result skipped: ", partner, " / ", timepoint)
        next
      }
      
      result <- tryCatch(
        run_partner_timepoint(
          partner = partner,
          signal_matrix = signal_matrix,
          doses = doses,
          output_prefix = output_prefix
        ),
        error = function(e) {
          error_file <- paste0(output_prefix, "_ERROR.txt")
          writeLines(conditionMessage(e), error_file)
          warning(
            partner, " / ", permutation_id, " / ", timepoint,
            " failed: ", conditionMessage(e)
          )
          NULL
        }
      )
      
      if (is.null(result)) {
        summary_rows[[summary_index]] <- data.frame(
          permutation = permutation_id,
          partner = partner,
          timepoint = timepoint,
          normalization = partner_objects[[partner]]$method,
          n_samples = length(idx),
          n_unique_doses = length(unique(doses)),
          n_selected = NA_integer_,
          status = "ERROR",
          stringsAsFactors = FALSE
        )
        summary_index <- summary_index + 1L
        next
      }
      
      run_info <- data.frame(
        permutation = permutation_id,
        partner = partner,
        timepoint = timepoint,
        normalization = partner_objects[[partner]]$method,
        FDR = PREFILTER_FDR,
        n_samples = length(idx),
        n_unique_doses = length(unique(doses)),
        stringsAsFactors = FALSE
      )
      
      save_selection(
        result$selection,
        output_prefix,
        run_info
      )
      
      summary_rows[[summary_index]] <- data.frame(
        permutation = permutation_id,
        partner = partner,
        timepoint = timepoint,
        normalization = partner_objects[[partner]]$method,
        n_samples = length(idx),
        n_unique_doses = length(unique(doses)),
        n_selected = result$n_selected,
        status = "OK",
        stringsAsFactors = FALSE
      )
      summary_index <- summary_index + 1L
    }
  }
  
  # Continuously update the summary so progress survives an interrupted run.
  if (length(summary_rows) > 0L) {
    write.csv(
      do.call(rbind, summary_rows),
      file.path(OUTPUT_DIR, "quadratic_prefilter_summary.csv"),
      row.names = FALSE
    )
  }
  
  gc(verbose = FALSE)
}

message(
  "Finished. Summary written to: ",
  file.path(OUTPUT_DIR, "quadratic_prefilter_summary.csv")
)
