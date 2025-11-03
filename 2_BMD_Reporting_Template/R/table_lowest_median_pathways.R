
# Top N lowest-median-BMD pathways per timepoint
top_pathways_per_timepoint <- function(
    df,
    n = 5,
    digits = 3,
    caption_prefix = "Top pathways with lowest median BMD per timepoint",
    return_data = FALSE
) {
  stopifnot(is.data.frame(df), all(c("pathway", "timepoint", "bmd") %in% names(df)))
  
  # Prepare & clean
  df <- df %>%
    dplyr::select(dplyr::any_of(c("timepoint", "pathway", "bmd", "bmdl", "bmdu"))) %>%
    dplyr::filter(!is.na(bmd)) %>%
    dplyr::mutate(
      timepoint = as.character(timepoint),
      bmd = as.numeric(bmd)
    )
  
  # Order timepoints nicely
  tps_ord <- tryCatch(levels(order_tp(unique(df$timepoint))), error = function(e) unique(df$timepoint))
  
  # Get top N per timepoint
  top_tbl <- df %>%
    dplyr::group_by(timepoint) %>%
    dplyr::arrange(bmd, .by_group = TRUE) %>%
    dplyr::slice_head(n = n) %>%
    dplyr::ungroup()
  
  # Round numeric columns
  num_cols <- names(top_tbl)[sapply(top_tbl, is.numeric)]
  top_tbl[num_cols] <- lapply(top_tbl[num_cols], function(x) round(x, digits))
  
  if (return_data) return(top_tbl)
  
  # Pretty Rmd output: one table per timepoint
  purrr::walk(tps_ord, function(tp) {
    sub_df <- dplyr::filter(top_tbl, timepoint == tp)
    if (nrow(sub_df) == 0) return(NULL)
    print(
      knitr::kable(
        sub_df %>% dplyr::select(-timepoint),
        caption = glue::glue("{caption_prefix} — {tp}"),
        format = "html",
        align = c("l", rep("r", ncol(sub_df) - 2))
      ) |>
        kableExtra::kable_styling(
          bootstrap_options = c("striped", "hover", "condensed"),
          full_width = FALSE
        )
    )
  })
}
