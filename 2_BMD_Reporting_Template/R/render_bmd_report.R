

render_bmd_report <- function(
    rmd_path,
    reports_dir = "reports",
    output_basename = NULL,
    params_override = list(),
    quiet = TRUE
) {
  stopifnot(file.exists(rmd_path))
  
  rmd_path <- normalizePath(rmd_path, winslash = "/", mustWork = TRUE)
  rmd_dir  <- dirname(rmd_path)
  
  # --- unwrap params like list(value=...) recursively ---
  unwrap_param_values <- function(x) {
    if (!is.list(x)) return(x)
    if (!is.null(x$value) && ("value" %in% names(x))) return(unwrap_param_values(x$value))
    out <- lapply(x, unwrap_param_values)
    if (!is.null(names(x))) names(out) <- names(x)
    out
  }
  
  default_params_raw <- rmarkdown::yaml_front_matter(rmd_path)$params
  if (is.null(default_params_raw)) default_params_raw <- list()
  default_params <- unwrap_param_values(default_params_raw)
  
  params <- modifyList(default_params, params_override)
  
  if (is.null(output_basename) || !nzchar(output_basename)) {
    base <- tools::file_path_sans_ext(basename(rmd_path))
    study_id <- if (!is.null(params$study_id)) params$study_id else "study"
    stamp <- format(Sys.Date(), "%Y%m%d")
    output_basename <- paste0(base, "_", study_id, "_", stamp, ".html")
  }
  if (!grepl("\\.html?$", output_basename, ignore.case = TRUE)) {
    output_basename <- paste0(output_basename, ".html")
  }
  
  out_dir <- file.path(rmd_dir, reports_dir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Make relative paths behave like knitting from the Rmd location
  old_root <- knitr::opts_knit$get("root.dir")
  on.exit(knitr::opts_knit$set(root.dir = old_root), add = TRUE)
  knitr::opts_knit$set(root.dir = rmd_dir)
  
  rendered_path <- rmarkdown::render(
    input = rmd_path,
    output_format = "html_document",
    output_file = output_basename,
    output_dir = out_dir,
    params = params,
    envir = globalenv(),
    quiet = quiet
  )
  
  normalizePath(rendered_path, winslash = "/", mustWork = TRUE)
}
