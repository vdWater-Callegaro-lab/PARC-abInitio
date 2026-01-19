
# Top N lowest-pthway-BMD per timepoint

# ---- helper: compute top N per timepoint (shared by both functions) ----
compute_top_pathways <- function(df, n = 5, digits = 3) { 
  stopifnot(is.data.frame(df), all(c("pathway_id", "timepoint", "bmd") %in% names(df)))
  
  df <- df %>%
    dplyr::select(dplyr::any_of(c("timepoint", "pathway_id", "bmd", "bmdl", "bmdu"))) %>%
    dplyr::filter(!is.na(bmd))
  
  if (!"bmdl" %in% names(df)) df$bmdl <- NA_real_
  if (!"bmdu" %in% names(df)) df$bmdu <- NA_real_
  
  top_tbl <- df %>%
    dplyr::mutate(
      timepoint  = as.character(timepoint),
      pathway_id = as.character(pathway_id),
      bmd  = as.numeric(bmd),
      bmdl = as.numeric(bmdl),
      bmdu = as.numeric(bmdu)
    ) %>%
    dplyr::group_by(timepoint) %>%
    dplyr::arrange(bmd, .by_group = TRUE) %>%
    dplyr::slice_head(n = n) %>%
    dplyr::ungroup()
  
  # round numeric columns for display / export
  num_cols <- names(top_tbl)[vapply(top_tbl, is.numeric, logical(1))]
  top_tbl[num_cols] <- lapply(top_tbl[num_cols], function(x) round(x, digits))
  
  top_tbl %>%
    dplyr::select(timepoint, pathway_id, dplyr::any_of(c("bmd", "bmdl", "bmdu")))
}


# ---- 1) render interactive tables (one per timepoint) ----
render_top_pathway_tables <- function( 
    df,
    n = 5,
    digits = 3,
    caption_prefix = "Top pathways with lowest BMD per timepoint",
    page_length = 15,
    return_tables = FALSE
) {
  stopifnot(is.data.frame(df), all(c("pathway_id", "timepoint", "bmd") %in% names(df)))
  
  df <- df %>%
    dplyr::select(dplyr::any_of(c("timepoint", "pathway_id", "bmd", "bmdl", "bmdu"))) %>%
    dplyr::filter(!is.na(bmd))
  
  if (!"bmdl" %in% names(df)) df$bmdl <- NA_real_
  if (!"bmdu" %in% names(df)) df$bmdu <- NA_real_
  
  out <- df %>%
    dplyr::mutate(
      timepoint  = as.character(timepoint),
      pathway_id = as.character(pathway_id),
      bmd  = as.numeric(bmd),
      bmdl = as.numeric(bmdl),
      bmdu = as.numeric(bmdu)
    ) %>%
    dplyr::group_by(timepoint) %>%
    dplyr::arrange(bmd, .by_group = TRUE) %>%
    dplyr::slice_head(n = n) %>%
    dplyr::ungroup()
  
  num_cols <- names(out)[vapply(out, is.numeric, logical(1))]
  out[num_cols] <- lapply(out[num_cols], function(x) round(x, digits))
  
  out <- out %>%
    dplyr::select(timepoint, pathway_id, dplyr::any_of(c("bmd", "bmdl", "bmdu")))
  
  # ---- ORDER timepoints numerically (8h, 16h, 24h, 48h, ...) ----
  tp_raw <- unique(out$timepoint)
  
  # Extract first number from each timepoint string (works for "8h", "8hr", "TP8", etc.)
  tp_num <- suppressWarnings(as.numeric(gsub(".*?(\\d+).*", "\\1", tp_raw)))
  
  # If extraction fails for some values, keep them at the end in original order
  tp_ordered <- c(tp_raw[order(tp_num, na.last = TRUE)])
  
  # ---- Make one DT table per timepoint ----
  tables <- lapply(tp_ordered, function(tp) {
    sub <- dplyr::filter(out, timepoint == tp) %>%
      dplyr::select(-timepoint) %>%
      dplyr::rename(Pathway = pathway_id, BMD = bmd, BMDL = bmdl, BMDU = bmdu)
    
    DT::datatable(
      sub,
      rownames = FALSE,
      options = list(
        pageLength = page_length,
        lengthMenu = c(15, 25, 50, 100),
        searching = TRUE,
        ordering = TRUE
      ),
      caption = htmltools::tags$caption(
        style = "caption-side: top; text-align: left; font-weight: bold;",
        paste0(caption_prefix, " — ", tp)
      )
    )
  })
  
  names(tables) <- tp_ordered
  if (return_tables) return(tables)
  
  # Return as a single HTML object so it renders correctly in Rmd
  htmltools::tagList(
    lapply(tp_ordered, function(tp) {
      htmltools::tagList(
        htmltools::tags$h3(tp),
        tables[[tp]],
        htmltools::tags$br()
      )
    })
  )
}




# ---- 2) save combined overview table (all timepoints) ----
save_top_pathway_overview <- function( 
    df,
    output_dir,
    filename = "top_pathways_overview.txt",
    n = 5,
    digits = 3
) {
  stopifnot(is.data.frame(df), all(c("pathway_id", "timepoint", "bmd") %in% names(df)))
  
  # keep relevant cols
  df <- df %>%
    dplyr::select(dplyr::any_of(c("timepoint", "pathway_id", "bmd", "bmdl", "bmdu"))) %>%
    dplyr::filter(!is.na(bmd))
  
  # ensure optional cols exist
  if (!"bmdl" %in% names(df)) df$bmdl <- NA_real_
  if (!"bmdu" %in% names(df)) df$bmdu <- NA_real_
  
  # top N per timepoint
  overview_tbl <- df %>%
    dplyr::mutate(
      timepoint  = as.character(timepoint),
      pathway_id = as.character(pathway_id),
      bmd  = as.numeric(bmd),
      bmdl = as.numeric(bmdl),
      bmdu = as.numeric(bmdu)
    ) %>%
    dplyr::group_by(timepoint) %>%
    dplyr::arrange(bmd, .by_group = TRUE) %>%
    dplyr::slice_head(n = n) %>%
    dplyr::ungroup()
  
  # round numeric columns
  num_cols <- names(overview_tbl)[vapply(overview_tbl, is.numeric, logical(1))]
  overview_tbl[num_cols] <- lapply(overview_tbl[num_cols], function(x) round(x, digits))
  
  # final overview (keep order)
  overview_tbl <- overview_tbl %>%
    dplyr::select(timepoint, pathway_id, dplyr::any_of(c("bmd", "bmdl", "bmdu")))
  
  # save
  out_path <- file.path(output_dir, filename)
  data.table::fwrite(
    overview_tbl,
    file = out_path,
    sep = "\t"
  )
  
  message("Saved overview table to: ", out_path)
  invisible(overview_tbl)
}



make_top_pw_pertp <- function(df, k = 4, digits = 3) {
  stopifnot(is.data.frame(df), all(c("pathway_id", "timepoint", "bmd") %in% names(df)))
  stopifnot(is.numeric(k), length(k) == 1, k >= 1)
  
  df <- df %>%
    dplyr::select(dplyr::any_of(c("timepoint", "pathway_id", "bmd", "bmdl", "bmdu"))) %>%
    dplyr::filter(!is.na(bmd))
  
  if (!"bmdl" %in% names(df)) df$bmdl <- NA_real_
  if (!"bmdu" %in% names(df)) df$bmdu <- NA_real_
  
  top_pw_pertp <- df %>%
    dplyr::mutate(
      timepoint  = as.character(timepoint),
      pathway_id = as.character(pathway_id),
      bmd  = as.numeric(bmd),
      bmdl = as.numeric(bmdl),
      bmdu = as.numeric(bmdu)
    ) %>%
    dplyr::group_by(timepoint) %>%
    dplyr::arrange(bmd, .by_group = TRUE) %>%
    dplyr::slice_head(n = k) %>%
    dplyr::ungroup()
  
  # round numeric columns (for display/use)
  num_cols <- names(top_pw_pertp)[vapply(top_pw_pertp, is.numeric, logical(1))]
  top_pw_pertp[num_cols] <- lapply(top_pw_pertp[num_cols], function(x) round(x, digits))
  
  top_pw_pertp %>%
    dplyr::select(timepoint, pathway_id, dplyr::any_of(c("bmd", "bmdl", "bmdu")))
}






