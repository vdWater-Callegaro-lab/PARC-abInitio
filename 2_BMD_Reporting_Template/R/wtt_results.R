
# WTT - plot dose responsive counts per timepoint (standardized input)
doseresponsive_plot <- function(
    WTT,
    level = "genes",
    p_col = "p",                 # choose column for significance testing (either p or p_adj)
    p_thresh = 0.05,                   # threshold for significance 
    palette = c(upregulated = "firebrick", downregulated = "darkblue")
) {
  stopifnot(is.data.frame(WTT), all(c("timepoint", "max_fold_change") %in% names(WTT)))
  if (!p_col %in% names(WTT)) stop(glue::glue("Column '{p_col}' not found in WTT."))
  
  # Determine regulation (up or down)
  WTT <- WTT %>%
    mutate(
      regulation = ifelse(max_fold_change > 0, "upregulated", "downregulated")
    )
  
  # Filter significant results
  df_sig <- WTT %>% filter(.data[[p_col]] < p_thresh)
  
  # Ensure consistent timepoint ordering
  tps_all <- unique(as.character(WTT$timepoint))
  tps_ord <- tryCatch(levels(order_tp(tps_all)), error = function(e) tps_all)
  
  # Build full grid to include zeros
  regs <- c("upregulated", "downregulated")
  all_combos <- tidyr::expand_grid(timepoint = tps_ord, regulation = regs)
  
  counts <- df_sig %>%
    group_by(timepoint, regulation) %>%
    summarise(count = n(), .groups = "drop") %>%
    right_join(all_combos, by = c("timepoint", "regulation")) %>%
    mutate(
      count = replace_na(count, 0L),
      timepoint = factor(timepoint, levels = tps_ord)
    )
  
  ggplot(counts, aes(x = timepoint, y = count, fill = regulation)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = palette, name = "Regulation") +
    labs(
      title = paste0("Number of dose responsive ", level, " per timepoint"),
      subtitle = glue::glue("Significance: {p_col} < {p_thresh}"),
      x = NULL,
      y = paste0("Number of dose responsive ", level)
    ) +
    theme_minimal(base_size = 11)
}




# Dose-responsive counts per timepoint
dose_responsive_table <- function(
    WTT,
    p_col   = "p",      # or "p_adj"
    p_thresh = 0.05,
    timepoints = NULL,        # optional vector to force/include specific TPs (even with zero hits)
    caption = NULL,           # optional table caption
    return_data = FALSE       # set TRUE to get the underlying data.frame instead of a kable
) {
  stopifnot(is.data.frame(WTT), all(c("timepoint", "max_fold_change") %in% names(WTT)))
  if (!p_col %in% names(WTT)) stop(glue::glue("Column '{p_col}' not found in WTT."))
  
  # determine regulation and filter by significance
  df <- WTT %>%
    dplyr::mutate(regulation = ifelse(max_fold_change > 0, "upregulated", "downregulated")) %>%
    dplyr::filter(.data[[p_col]] < p_thresh)
  
  # timepoints to display (order naturally if possible)
  tps_all <- if (is.null(timepoints)) unique(as.character(WTT$timepoint)) else as.character(timepoints)
  tps_ord <- tryCatch(levels(order_tp(tps_all)), error = function(e) tps_all)
  
  # full grid to include zeros
  regs <- c("upregulated", "downregulated")
  all_grid <- tidyr::expand_grid(timepoint = tps_ord, regulation = regs)
  
  # counts + totals
  counts <- df %>%
    dplyr::count(timepoint, regulation, name = "n") %>%
    dplyr::right_join(all_grid, by = c("timepoint", "regulation")) %>%
    dplyr::mutate(n = tidyr::replace_na(n, 0L),
                  timepoint = factor(timepoint, levels = tps_ord)) %>%
    tidyr::pivot_wider(names_from = regulation, values_from = n, values_fill = 0) %>%
    dplyr::arrange(timepoint) %>%
    dplyr::mutate(Total = upregulated + downregulated) %>%
    dplyr::rename(`Upregulated` = upregulated, `Downregulated` = downregulated) %>%
    dplyr::mutate(timepoint = as.character(timepoint))
  
  if (return_data) return(counts)
  
  # pretty Rmd table
  knitr::kable(
    counts,
    caption = if (!is.null(caption)) caption else
      glue::glue("Dose-responsive {ifelse('gene' %in% names(WTT), 'genes', 'entities')} per timepoint ({p_col} < {p_thresh})"),
    align = c("l", "r", "r", "r"),
    format = "html"
  ) |>
    kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE) |>
    kableExtra::add_header_above(c(" " = 1, "Dose-responsive counts" = 3))
}
