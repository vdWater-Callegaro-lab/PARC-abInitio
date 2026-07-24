#!/usr/bin/env Rscript

# ==============================================================================
# Sciensano empirical-null preprocessing and BMDExpress upload generation
#
# For every permuted metadata file:
#   1. Align unchanged raw counts with the permuted metadata.
#   2. Reproduce Sciensano's RODAF-style CPM prefilter:
#        - calculate CPM from raw library sizes
#        - define conditions by permuted dose and timepoint
#        - retain a probe when CPM >= 1 in at least 75% of samples in
#          at least one condition
#   3. Join the manifest and sum probes by the Sciensano gene identifier.
#   4. Apply DESeq2 median-ratio normalization.
#   5. Replace normalized zeros by 1 and log2-transform.
#   6. Write one tab-delimited BMDExpress upload file for each
#      permutation x timepoint.
#
# No samples are removed.
# ==============================================================================


# ---- Package checks -----------------------------------------------------------

required_packages <- c("DESeq2")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install the following package(s) before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(DESeq2)
})


# ---- USER SETTINGS ------------------------------------------------------------

COUNT_FILE <- "data/EUT046/EUT046_counts_per_gene_per_sample_raw_final.csv"

MANIFEST_FILE <- "data/Manifests/Temposeq_manifest_Human_WT_1.2_realignment_96percent_2023-05-15.txt"

PERMUTED_METADATA_DIR <- "output/EUT046/permuted_metadata"

OUTPUT_DIR <- "output/EUT046/Sciensano_BMDExpress_uploads"

METADATA_PATTERN <- "\\.(csv|tsv|txt)$"

# Candidate metadata columns. The first matching name is used.
# Candidate metadata columns. The first matching name is used.
SAMPLE_COLUMN_CANDIDATES <- "SAMPLE_ID"

DOSE_COLUMN_CANDIDATES <- "CONCENTRATION"

TIME_COLUMN_CANDIDATES <- "TIME"

REPLICATE_COLUMN_CANDIDATES <- "REPLICATE"

COUNT_ID_COLUMN_CANDIDATES <- "Column1"

MANIFEST_PROBE_COLUMN_CANDIDATES <- c(
  "Probe_ID", "PROBE_ID", "probe_id", "ProbeID"
)

# The Sciensano BMDExpress files displayed identifiers such as A2M_2.
# Therefore gene_symbol_entrez_id is preferred over gene_symbol.
MANIFEST_GENE_COLUMN_CANDIDATES <- c(
  "gene_symbol_entrez_id",
  "gene_symbol",
  "Gene_Symbol_TS",
  "entrez_id"
)

CPM_THRESHOLD <- 1
MIN_FRACTION_WITHIN_CONDITION <- 0.75

# Sciensano's notebook replaces normalized zeros with 1 before log2.
ZERO_REPLACEMENT <- 1

# BMDExpress first-column header used in their upload notebook.
BMD_FIRST_COLUMN_HEADER <- "0_0_Dose"

# Number formatting for expression values and doses.
SIGNIFICANT_DIGITS <- 10

OVERWRITE <- FALSE


# ---- General helpers ----------------------------------------------------------

safe_dir_create <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}


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
  
  original <- trimws(as.character(x))
  original[original %in% c("", "NA", "NaN", "null", "NULL")] <- NA_character_
  
  cleaned <- gsub(",", ".", original, fixed = TRUE)
  cleaned <- gsub("[^0-9eE+.\\-]", "", cleaned, perl = TRUE)
  
  out <- suppressWarnings(as.numeric(cleaned))
  
  if (anyNA(out)) {
    bad <- unique(original[is.na(out)])
    stop(
      "Dose conversion failed for: ",
      paste(bad, collapse = ", ")
    )
  }
  
  out
}


standardize_time <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("(?i)hours?$", "h", x, perl = TRUE)
  x <- sub("(?i)hrs?$", "h", x, perl = TRUE)
  
  # Numeric timepoints are given an h suffix only for filenames/condition labels.
  numeric_only <- grepl("^[0-9.]+$", x)
  x[numeric_only] <- paste0(x[numeric_only], "h")
  
  x
}


safe_filename <- function(x) {
  gsub("[^A-Za-z0-9._-]", "_", x)
}


format_numeric <- function(x) {
  format(
    signif(as.numeric(x), SIGNIFICANT_DIGITS),
    scientific = FALSE,
    trim = TRUE
  )
}


# ---- Data import ---------------------------------------------------------------

read_count_matrix <- function(path) {
  if (!file.exists(path)) {
    stop("COUNT_FILE does not exist: ", path)
  }
  
  x <- read_table_auto(path, check.names = FALSE)
  
  id_col <- first_existing_name(
    names(x),
    COUNT_ID_COLUMN_CANDIDATES,
    "the raw-count probe identifier column"
  )
  
  ids <- as.character(x[[id_col]])
  x[[id_col]] <- NULL
  
  x[] <- lapply(x, function(z) {
    out <- suppressWarnings(as.numeric(as.character(z)))
    
    if (anyNA(out) && !all(is.na(z))) {
      stop("Non-numeric values were found in a raw-count sample column.")
    }
    
    out
  })
  
  mat <- as.matrix(x)
  rownames(mat) <- ids
  
  if (anyDuplicated(rownames(mat))) {
    stop(
      "Raw count file contains duplicated probe IDs. ",
      "Sciensano's filtering expects one row per probe."
    )
  }
  
  if (anyDuplicated(colnames(mat))) {
    stop("Raw count file contains duplicated sample names.")
  }
  
  if (any(mat < 0, na.rm = TRUE)) {
    stop("Raw count matrix contains negative values.")
  }
  
  mat
}


read_manifest_map <- function(path) {
  if (!file.exists(path)) {
    stop("MANIFEST_FILE does not exist: ", path)
  }
  
  manifest <- read_table_auto(path, check.names = FALSE)
  
  probe_col <- first_existing_name(
    names(manifest),
    MANIFEST_PROBE_COLUMN_CANDIDATES,
    "the manifest probe-ID column"
  )
  
  gene_col <- first_existing_name(
    names(manifest),
    MANIFEST_GENE_COLUMN_CANDIDATES,
    "the manifest gene identifier column"
  )
  
  map <- data.frame(
    probe_id = as.character(manifest[[probe_col]]),
    gene_id = as.character(manifest[[gene_col]]),
    stringsAsFactors = FALSE
  )
  
  map$gene_id[
    map$gene_id %in% c("", "NA", "NaN", "None", "null", "NULL")
  ] <- NA_character_
  
  # Python/pandas may render integer-like Entrez values as "2.0".
  map$gene_id <- sub("\\.0$", "", map$gene_id)
  
  map
}


infer_replicate <- function(sample_names) {
  out <- sub(
    ".*(?:_Rep|_REP|_rep|_R|_r)([0-9]+)$",
    "\\1",
    sample_names,
    perl = TRUE
  )
  
  failed <- identical(out, sample_names)
  
  if (failed) {
    return(seq_along(sample_names))
  }
  
  out
}


read_permuted_metadata <- function(path) {
  meta <- read_table_auto(path, check.names = FALSE)
  
  sample_col <- first_existing_name(
    names(meta),
    SAMPLE_COLUMN_CANDIDATES,
    "the metadata sample column"
  )
  
  dose_col <- first_existing_name(
    names(meta),
    DOSE_COLUMN_CANDIDATES,
    "the metadata dose/concentration column"
  )
  
  time_col <- first_existing_name(
    names(meta),
    TIME_COLUMN_CANDIDATES,
    "the metadata timepoint column"
  )
  
  replicate_hits <- REPLICATE_COLUMN_CANDIDATES[
    REPLICATE_COLUMN_CANDIDATES %in% names(meta)
  ]
  
  sample_names <- as.character(meta[[sample_col]])
  
  replicate_values <- if (length(replicate_hits) > 0L) {
    as.character(meta[[replicate_hits[[1L]]]])
  } else {
    infer_replicate(sample_names)
  }
  
  out <- data.frame(
    sample = sample_names,
    dose = clean_numeric_dose(meta[[dose_col]]),
    timepoint = standardize_time(meta[[time_col]]),
    replicate = replicate_values,
    stringsAsFactors = FALSE
  )
  
  if (anyDuplicated(out$sample)) {
    stop("Duplicated samples in metadata: ", basename(path))
  }
  
  if (anyNA(out$sample) || any(out$sample == "")) {
    stop("Missing sample names in metadata: ", basename(path))
  }
  
  if (anyNA(out$dose)) {
    stop("Missing dose values in metadata: ", basename(path))
  }
  
  if (anyNA(out$timepoint) || any(out$timepoint == "")) {
    stop("Missing timepoints in metadata: ", basename(path))
  }
  
  out
}


align_counts_and_metadata <- function(counts, metadata) {
  missing_in_counts <- setdiff(metadata$sample, colnames(counts))
  extra_in_counts <- setdiff(colnames(counts), metadata$sample)
  
  if (length(missing_in_counts) > 0L) {
    stop(
      "Metadata samples absent from count matrix: ",
      paste(missing_in_counts, collapse = ", ")
    )
  }
  
  # User confirmed that Sciensano did not remove samples.
  if (length(extra_in_counts) > 0L) {
    stop(
      "Count-matrix samples absent from metadata: ",
      paste(extra_in_counts, collapse = ", ")
    )
  }
  
  counts <- counts[, metadata$sample, drop = FALSE]
  
  stopifnot(identical(colnames(counts), metadata$sample))
  
  list(counts = counts, metadata = metadata)
}


# ---- Sciensano preprocessing --------------------------------------------------

calculate_cpm <- function(counts) {
  library_sizes <- colSums(counts)
  
  if (any(library_sizes <= 0)) {
    stop(
      "One or more samples have a non-positive raw library size: ",
      paste(colnames(counts)[library_sizes <= 0], collapse = ", ")
    )
  }
  
  sweep(counts, 2L, library_sizes, "/") * 1e6
}


rodaf_keep_probes <- function(counts, metadata) {
  cpm <- calculate_cpm(counts)
  
  # This reproduces the Python condition construction:
  # samples are grouped by the sample-name prefix before "_Rep".
  # For permuted analyses, the equivalent condition is permuted dose + timepoint.
  condition <- interaction(
    metadata$timepoint,
    metadata$dose,
    drop = TRUE,
    lex.order = TRUE
  )
  
  condition_levels <- levels(condition)
  
  pass_by_condition <- vapply(
    condition_levels,
    function(condition_level) {
      sample_idx <- which(condition == condition_level)
      required_n <- ceiling(
        MIN_FRACTION_WITHIN_CONDITION * length(sample_idx)
      )
      
      rowSums(
        cpm[, sample_idx, drop = FALSE] >= CPM_THRESHOLD
      ) >= required_n
    },
    logical(nrow(cpm))
  )
  
  if (is.null(dim(pass_by_condition))) {
    pass_by_condition <- matrix(
      pass_by_condition,
      ncol = 1L,
      dimnames = list(
        rownames(cpm),
        condition_levels
      )
    )
  }
  
  keep <- rowSums(pass_by_condition) >= 1L
  
  list(
    keep = keep,
    pass_by_condition = pass_by_condition,
    cpm = cpm
  )
}


aggregate_probes_to_genes <- function(filtered_counts, manifest_map) {
  gene_ids <- manifest_map$gene_id[
    match(rownames(filtered_counts), manifest_map$probe_id)
  ]
  
  keep_mapped <- !is.na(gene_ids) & gene_ids != ""
  
  if (!any(keep_mapped)) {
    stop(
      "No retained probes could be mapped to a Sciensano gene identifier."
    )
  }
  
  mapped_counts <- filtered_counts[keep_mapped, , drop = FALSE]
  mapped_gene_ids <- gene_ids[keep_mapped]
  
  summed <- rowsum(
    mapped_counts,
    group = mapped_gene_ids,
    reorder = FALSE,
    na.rm = TRUE
  )
  
  storage.mode(summed) <- "integer"
  summed
}


deseq2_normalized_counts <- function(summed_counts, metadata) {
  # Sciensano used DESeq2 and then counts(dds, normalized = TRUE).
  #
  # Size-factor estimation does not depend on the design coefficients.
  # A design of ~1 therefore gives the same median-ratio normalized counts
  # while avoiding failures caused by label permutations creating aliased
  # replicate/treatment combinations.
  col_data <- data.frame(
    row.names = metadata$sample,
    replicate = factor(metadata$replicate),
    timepoint = factor(metadata$timepoint),
    dose = factor(metadata$dose)
  )
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(summed_counts),
    colData = col_data,
    design = ~ 1
  )
  
  dds <- DESeq2::estimateSizeFactors(dds)
  
  normalized <- DESeq2::counts(
    dds,
    normalized = TRUE
  )
  
  normalized
}


log2_for_bmdexpress <- function(normalized_counts) {
  # Exact behavior in the Sciensano formatting notebook:
  # replace 0 by 1, then take log2.
  normalized_counts[normalized_counts == 0] <- ZERO_REPLACEMENT
  log2(normalized_counts)
}


# ---- BMDExpress formatting ----------------------------------------------------

make_bmdexpress_table <- function(
    log2_matrix,
    metadata_timepoint
) {
  sample_names <- metadata_timepoint$sample
  doses <- metadata_timepoint$dose
  
  expr <- log2_matrix[, sample_names, drop = FALSE]
  
  # BMDExpress layout used by Sciensano:
  #   first row:  Dose | dose per sample
  #   later rows: gene | expression per sample
  dose_row <- data.frame(
    feature = "Dose",
    t(format_numeric(doses)),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  expression_rows <- data.frame(
    feature = rownames(expr),
    expr,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  names(dose_row) <- c(
    BMD_FIRST_COLUMN_HEADER,
    sample_names
  )
  
  names(expression_rows) <- c(
    BMD_FIRST_COLUMN_HEADER,
    sample_names
  )
  
  out <- rbind(
    dose_row,
    expression_rows
  )
  
  out
}


validate_bmdexpress_table <- function(x, expected_samples) {
  if (!identical(names(x)[-1L], expected_samples)) {
    stop("BMDExpress table sample order is incorrect.")
  }
  
  if (x[[1L]][1L] != "Dose") {
    stop("The first BMDExpress data row is not the Dose row.")
  }
  
  if (anyDuplicated(x[[1L]][-1L])) {
    warning(
      "Duplicated feature identifiers are present in the BMDExpress table."
    )
  }
  
  invisible(TRUE)
}


write_bmdexpress_file <- function(x, path) {
  write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  )
}


# ---- Main ---------------------------------------------------------------------

safe_dir_create(OUTPUT_DIR)

metadata_files <- list.files(
  PERMUTED_METADATA_DIR,
  pattern = METADATA_PATTERN,
  full.names = TRUE,
  ignore.case = TRUE
)

metadata_files <- sort(metadata_files)

if (length(metadata_files) == 0L) {
  stop(
    "No permuted metadata files found in: ",
    PERMUTED_METADATA_DIR
  )
}

if (length(metadata_files) != 100L) {
  warning(
    "Expected 100 permuted metadata files but found ",
    length(metadata_files),
    ". All files found will be processed."
  )
}

raw_counts <- read_count_matrix(COUNT_FILE)
manifest_map <- read_manifest_map(MANIFEST_FILE)

summary_rows <- list()
summary_i <- 1L

for (metadata_file in metadata_files) {
  permutation_id <- tools::file_path_sans_ext(
    basename(metadata_file)
  )
  
  message("Permutation: ", permutation_id)
  
  audit_dir <- file.path(
    OUTPUT_DIR,
    "_preprocessing_audit"
  )
  safe_dir_create(audit_dir)
  
  metadata <- read_permuted_metadata(metadata_file)
  
  aligned <- align_counts_and_metadata(
    raw_counts,
    metadata
  )
  
  counts_aligned <- aligned$counts
  metadata <- aligned$metadata
  
  error_path <- file.path(
    audit_dir,
    paste0(permutation_id, "_ERROR.txt")
  )
  
  permutation_result <- tryCatch(
    {
      filter_result <- rodaf_keep_probes(
        counts_aligned,
        metadata
      )
      
      filtered_counts <- counts_aligned[
        filter_result$keep,
        ,
        drop = FALSE
      ]
      
      if (nrow(filtered_counts) == 0L) {
        stop("RODAF-style filtering retained zero probes.")
      }
      
      summed_counts <- aggregate_probes_to_genes(
        filtered_counts,
        manifest_map
      )
      
      normalized_counts <- deseq2_normalized_counts(
        summed_counts,
        metadata
      )
      
      log2_matrix <- log2_for_bmdexpress(
        normalized_counts
      )
      
      # Save intermediate matrices once per permutation for auditability.
      saveRDS(
        list(
          filtered_probe_counts = filtered_counts,
          summed_gene_counts = summed_counts,
          normalized_gene_counts = normalized_counts,
          log2_normalized_gene_counts = log2_matrix,
          metadata = metadata,
          filter_pass_by_condition = filter_result$pass_by_condition
        ),
        file.path(
          audit_dir,
          paste0(permutation_id, "_Sciensano_preprocessing.rds")
        )
      )
      
      list(
        filtered_counts = filtered_counts,
        summed_counts = summed_counts,
        normalized_counts = normalized_counts,
        log2_matrix = log2_matrix
      )
    },
    error = function(e) {
      writeLines(conditionMessage(e), error_path)
      warning(
        permutation_id,
        " failed: ",
        conditionMessage(e)
      )
      NULL
    }
  )
  
  if (is.null(permutation_result)) {
    summary_rows[[summary_i]] <- data.frame(
      permutation = permutation_id,
      timepoint = NA_character_,
      n_samples = NA_integer_,
      n_input_probes = nrow(counts_aligned),
      n_filtered_probes = NA_integer_,
      n_summed_genes = NA_integer_,
      output_file = NA_character_,
      status = "ERROR",
      stringsAsFactors = FALSE
    )
    summary_i <- summary_i + 1L
    next
  }
  
  for (timepoint_value in unique(metadata$timepoint)) {
    idx <- which(metadata$timepoint == timepoint_value)
    metadata_timepoint <- metadata[idx, , drop = FALSE]
    
    safe_timepoint <- safe_filename(timepoint_value)
    
    timepoint_dir <- file.path(
      OUTPUT_DIR,
      safe_timepoint
    )
    safe_dir_create(timepoint_dir)
    
    output_path <- file.path(
      timepoint_dir,
      paste0(
        permutation_id,
        "_",
        safe_timepoint,
        "_Sciensano_for_BMDExpress.txt"
      )
    )
    
    if (file.exists(output_path) && !OVERWRITE) {
      message("  Existing file skipped: ", basename(output_path))
      
      summary_rows[[summary_i]] <- data.frame(
        permutation = permutation_id,
        timepoint = timepoint_value,
        n_samples = nrow(metadata_timepoint),
        n_input_probes = nrow(counts_aligned),
        n_filtered_probes = nrow(permutation_result$filtered_counts),
        n_summed_genes = nrow(permutation_result$summed_counts),
        output_file = output_path,
        status = "EXISTING",
        stringsAsFactors = FALSE
      )
      summary_i <- summary_i + 1L
      next
    }
    
    bmd_table <- make_bmdexpress_table(
      permutation_result$log2_matrix,
      metadata_timepoint
    )
    
    validate_bmdexpress_table(
      bmd_table,
      metadata_timepoint$sample
    )
    
    write_bmdexpress_file(
      bmd_table,
      output_path
    )
    
    message("  Wrote: ", output_path)
    
    summary_rows[[summary_i]] <- data.frame(
      permutation = permutation_id,
      timepoint = timepoint_value,
      n_samples = nrow(metadata_timepoint),
      n_input_probes = nrow(counts_aligned),
      n_filtered_probes = nrow(permutation_result$filtered_counts),
      n_summed_genes = nrow(permutation_result$summed_counts),
      output_file = output_path,
      status = "OK",
      stringsAsFactors = FALSE
    )
    summary_i <- summary_i + 1L
  }
  
  write.csv(
    do.call(rbind, summary_rows),
    file.path(
      OUTPUT_DIR,
      "Sciensano_BMDExpress_upload_summary.csv"
    ),
    row.names = FALSE
  )
  
  gc(verbose = FALSE)
}

message(
  "Finished. Summary written to: ",
  file.path(
    OUTPUT_DIR,
    "Sciensano_BMDExpress_upload_summary.csv"
  )
)
