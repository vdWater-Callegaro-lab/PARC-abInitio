

# Read anything (.csv, .tsv, .txt, .xlsx). If sheet is needed, pass like "file.xlsx:Sheet1"
read_any <- function(path_like) {
  stopifnot(!is.null(path_like), length(path_like) == 1)
  
  # Allow sheet syntax: "file.xlsx:Sheet1"
  if (stringr::str_detect(path_like, ":")) {
    fp <- stringr::str_split_fixed(path_like, ":", 2)[, 1]
    sh <- stringr::str_split_fixed(path_like, ":", 2)[, 2]
  } else {
    fp <- path_like
    sh <- NA_character_
  }
  
  # Build path robustly
  base_dir <- tryCatch(
    normalizePath(as.character(input_dir), winslash = "/", mustWork = FALSE),
    error = function(e) as.character(input_dir)
  )
  full <- tryCatch(
    normalizePath(file.path(base_dir, fp), winslash = "/", mustWork = FALSE),
    error = function(e) file.path(base_dir, fp)
  )
  
  if (!file.exists(full)) {
    stop(glue::glue(
      "Missing file: {full}\n  - working dir: {getwd()}\n  - input_dir: {base_dir}"
    ))
  }
  
  ext <- tolower(tools::file_ext(full))
  
  if (ext %in% c("xlsx")) {
    return(readxl::read_xlsx(full, sheet = if (!is.na(sh)) sh else 1))
  }
  
  # For all delimited text (csv/tsv/txt/others), let fread auto-detect the separator,
  # handle multiple spaces/tabs, quoted fields, missing columns, etc.
  data.table::fread(
    full,
    sep = "auto",
    data.table = FALSE,
    na.strings = c("NA", "", "NaN"),
    check.names = T
  )
}

# Try to rename columns using a mapping list of synonyms
rename_using_mapping <- function(df, mapping) {
  df <- janitor::clean_names(df)
  current <- names(df)
  for (std in names(mapping)) {
    # mapping[[std]] is a character vector of candidate names (already clean_names style)
    candidates <- janitor::make_clean_names(mapping[[std]])
    hit <- intersect(candidates, current)
    if (length(hit) >= 1) {
      names(df)[match(hit[1], names(df))] <- std
      current <- names(df)
    }
  }
  df
}

# Ensure required columns exist
require_cols <- function(df, cols, context="table") {
  missing <- setdiff(cols, names(df))
  if (length(missing)) stop(glue::glue("Missing required columns in {context}: {paste(missing, collapse=", ")}"))
  invisible(df)
}

# Order factors nicely
order_tp <- function(x) {
  # Try to parse numeric hours if present
  if (all(grepl("^[0-9]+", x))) {
    ord <- order(as.numeric(x))
    factor(x, levels = unique(x[ord]))
  } else { factor(x) }
}
