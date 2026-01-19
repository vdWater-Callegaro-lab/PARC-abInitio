
choose_p_col <- function(TT, p_col = NULL) {
  if (!is.null(p_col) && p_col %in% names(TT)) return(p_col)
  if ("adj_p" %in% names(TT)) return("adj_p")
  if ("p" %in% names(TT)) return("p")
  stop("No valid p-value column found in TT. Expected 'adj_p' or 'p'.")
}

tp_num <- function(x) suppressWarnings(as.numeric(sub(".*?(\\d+).*", "\\1", as.character(x))))

order_timepoints <- function(tps) {
  tps <- as.character(tps)
  nums <- tp_num(tps)
  # numeric ones first (ordered), then any non-numeric labels in original order
  c(tps[order(nums, na.last = TRUE)])
}



# TT - plot dose responsive counts per timepoint (standardized input)
doseresponsive_plot <- function(
    TT,
    level = "genes",
    p_col = NULL,                 # NULL = auto choose adj_p if present else p
    p_thresh = 0.05,
    palette = c(upregulated = "firebrick", downregulated = "darkblue"),
    total_fill = "grey60"
) {
  stopifnot(is.data.frame(TT), "timepoint" %in% names(TT))
  p_col <- choose_p_col(TT, p_col)
  
  # significant rows
  df_sig <- TT %>% dplyr::filter(.data[[p_col]] < p_thresh)
  
  # order timepoints numerically
  tps_all <- unique(as.character(TT$timepoint))
  tps_ord <- order_timepoints(tps_all)
  
  # check if we can do up/down
  has_fc <- "max_fold_change" %in% names(TT) &&
    any(!is.na(suppressWarnings(as.numeric(TT$max_fold_change))))
  
  # ---- total only (no max_fold_change) ----
  if (!has_fc) {
    counts <- df_sig %>%
      dplyr::count(timepoint, name = "count") %>%
      dplyr::right_join(tibble::tibble(timepoint = tps_ord), by = "timepoint") %>%
      dplyr::mutate(
        count = tidyr::replace_na(count, 0L),
        timepoint = factor(timepoint, levels = tps_ord)
      )
    
    return(
      ggplot2::ggplot(counts, ggplot2::aes(x = timepoint, y = count)) +
        ggplot2::geom_col(fill = total_fill) +
        ggplot2::labs(
          title = paste0("Number of dose responsive ", level, " per timepoint"),
          subtitle = glue::glue("Significance: {p_col} < {p_thresh}"),
          x = NULL,
          y = paste0("Number of dose responsive ", level)
        ) +
        ggplot2::theme_minimal(base_size = 11)
    )
  }
  
  # ---- up/down breakdown (robust to existing TT$regulation) ----
  df_sig <- df_sig %>%
    dplyr::mutate(max_fold_change = as.numeric(max_fold_change)) %>%
    dplyr::select(-dplyr::any_of("regulation")) %>%
    dplyr::mutate(
      regulation = dplyr::if_else(
        !is.na(max_fold_change) & max_fold_change > 0,
        "upregulated",
        "downregulated"
      ),
      regulation = as.character(regulation)
    )
  
  regs <- c("upregulated", "downregulated")
  all_combos <- tidyr::expand_grid(timepoint = tps_ord, regulation = regs) %>%
    dplyr::mutate(
      timepoint = factor(timepoint, levels = tps_ord),
      regulation = as.character(regulation)
    )
  
  counts <- df_sig %>%
    dplyr::count(timepoint, regulation, name = "count") %>%
    dplyr::right_join(all_combos, by = c("timepoint", "regulation")) %>%
    dplyr::mutate(
      count = tidyr::replace_na(count, 0L),
      timepoint = factor(as.character(timepoint), levels = tps_ord)
    )
  
  ggplot2::ggplot(counts, ggplot2::aes(x = timepoint, y = count, fill = regulation)) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::scale_fill_manual(values = palette, name = "Regulation") +
    ggplot2::labs(
      title = paste0("Number of dose responsive ", level, " per timepoint"),
      subtitle = glue::glue("Significance: {p_col} < {p_thresh}"),
      x = NULL,
      y = paste0("Number of dose responsive ", level)
    ) +
    ggplot2::theme_minimal(base_size = 11)
}


# doseresponsive_plot <- function(
#     TT,
#     level = "genes",
#     p_col = "p",                 # choose column for significance testing (either p or p_adj)
#     p_thresh = 0.05,                   # threshold for significance 
#     palette = c(upregulated = "firebrick", downregulated = "darkblue")
# ) {
#   stopifnot(is.data.frame(TT), all(c("timepoint", "max_fold_change") %in% names(TT)))
#   if (!p_col %in% names(TT)) stop(glue::glue("Column '{p_col}' not found in TT."))
#   
#   # Determine regulation (up or down)
#   TT <- TT %>%
#     mutate(
#       regulation = ifelse(max_fold_change > 0, "upregulated", "downregulated")
#     )
#   
#   # Filter significant results
#   df_sig <- TT %>% filter(.data[[p_col]] < p_thresh)
#   
#   # Ensure consistent timepoint ordering
#   tps_all <- unique(as.character(TT$timepoint))
#   tps_ord <- tryCatch(levels(order_tp(tps_all)), error = function(e) tps_all)
#   
#   # Build full grid to include zeros
#   regs <- c("upregulated", "downregulated")
#   all_combos <- tidyr::expand_grid(timepoint = tps_ord, regulation = regs)
#   
#   counts <- df_sig %>%
#     group_by(timepoint, regulation) %>%
#     summarise(count = n(), .groups = "drop") %>%
#     right_join(all_combos, by = c("timepoint", "regulation")) %>%
#     mutate(
#       count = replace_na(count, 0L),
#       timepoint = factor(timepoint, levels = tps_ord)
#     )
#   
#   ggplot(counts, aes(x = timepoint, y = count, fill = regulation)) +
#     geom_col(position = "stack") +
#     scale_fill_manual(values = palette, name = "Regulation") +
#     labs(
#       title = paste0("Number of dose responsive ", level, " per timepoint"),
#       subtitle = glue::glue("Significance: {p_col} < {p_thresh}"),
#       x = NULL,
#       y = paste0("Number of dose responsive ", level)
#     ) +
#     theme_minimal(base_size = 11)
# }





# Dose-responsive counts per timepoint (simple text table)

dose_responsive_table <- function(
    TT,
    p_col = NULL,      # NULL = auto choose adj_p if present else p
    p_thresh = 0.05,
    timepoints = NULL, # optional: force/include specific TPs (even with zero hits)
    caption = NULL,
    return_data = FALSE
) {
  stopifnot(is.data.frame(TT), "timepoint" %in% names(TT))
  p_col <- choose_p_col(TT, p_col)
  
  df_sig <- TT %>% dplyr::filter(.data[[p_col]] < p_thresh)
  
  tps_all <- if (is.null(timepoints)) unique(as.character(TT$timepoint)) else as.character(timepoints)
  tps_ord <- order_timepoints(tps_all)
  
  has_fc <- "max_fold_change" %in% names(TT) &&
    any(!is.na(suppressWarnings(as.numeric(TT$max_fold_change))))
  
  # ---- total only ----
  if (!has_fc) {
    counts <- df_sig %>%
      dplyr::count(timepoint, name = "Total") %>%
      dplyr::right_join(tibble::tibble(timepoint = tps_ord), by = "timepoint") %>%
      dplyr::mutate(
        Total = tidyr::replace_na(Total, 0L),
        timepoint = factor(timepoint, levels = tps_ord)
      ) %>%
      dplyr::arrange(timepoint) %>%
      dplyr::mutate(timepoint = as.character(timepoint)) %>%
      as.data.frame(stringsAsFactors = FALSE)
    
    if (return_data) return(counts)
    
    if (is.null(caption)) {
      ent <- if ("gene" %in% names(TT)) "genes" else "entities"
      caption <- sprintf("Dose-responsive %s per timepoint (%s < %s)", ent, p_col, format(p_thresh))
    }
    
    names(counts)[names(counts) == "timepoint"] <- "Timepoint"
    return(knitr::kable(counts, caption = caption, align = "l"))
  }
  
  # ---- up/down breakdown (robust to existing TT$regulation) ----
  df <- df_sig %>%
    dplyr::mutate(max_fold_change = as.numeric(max_fold_change)) %>%
    dplyr::select(-dplyr::any_of("regulation")) %>%
    dplyr::mutate(
      regulation = dplyr::if_else(
        !is.na(max_fold_change) & max_fold_change > 0,
        "upregulated",
        "downregulated"
      ),
      regulation = as.character(regulation)
    )
  
  regs <- c("upregulated", "downregulated")
  all_grid <- tidyr::expand_grid(timepoint = tps_ord, regulation = regs) %>%
    dplyr::mutate(regulation = as.character(regulation))
  
  counts <- df %>%
    dplyr::count(timepoint, regulation, name = "n") %>%
    dplyr::right_join(all_grid, by = c("timepoint", "regulation")) %>%
    dplyr::mutate(
      n = tidyr::replace_na(n, 0L),
      timepoint = factor(as.character(timepoint), levels = tps_ord)
    ) %>%
    tidyr::pivot_wider(names_from = regulation, values_from = n, values_fill = 0) %>%
    dplyr::arrange(timepoint) %>%
    dplyr::mutate(Total = upregulated + downregulated) %>%
    dplyr::rename(Upregulated = upregulated, Downregulated = downregulated) %>%
    dplyr::mutate(timepoint = as.character(timepoint)) %>%
    as.data.frame(stringsAsFactors = FALSE)
  
  if (return_data) return(counts)
  
  if (is.null(caption)) {
    ent <- if ("gene" %in% names(TT)) "genes" else "entities"
    caption <- sprintf("Dose-responsive %s per timepoint (%s < %s)", ent, p_col, format(p_thresh))
  }
  
  names(counts)[names(counts) == "timepoint"] <- "Timepoint"
  knitr::kable(counts, caption = caption, align = "l")
}





# dose_responsive_table <- function(
#     TT,
#     p_col    = "p",      # or "p_adj"
#     p_thresh = 0.05,
#     timepoints = NULL,   # optional vector to force/include specific TPs (even with zero hits)
#     caption = NULL,      # optional caption printed above the table
#     return_data = FALSE  # TRUE returns the data.frame instead of printing
# ) {
#   stopifnot(is.data.frame(TT), all(c("timepoint", "max_fold_change") %in% names(TT)))
#   if (!p_col %in% names(TT)) stop(glue::glue("Column '{p_col}' not found in TT."))
#   
#   # determine regulation and filter by significance
#   df <- TT %>%
#     dplyr::mutate(regulation = ifelse(max_fold_change > 0, "upregulated", "downregulated")) %>%
#     dplyr::filter(.data[[p_col]] < p_thresh)
#   
#   # timepoints to display (order naturally if possible)
#   tps_all <- if (is.null(timepoints)) unique(as.character(TT$timepoint)) else as.character(timepoints)
#   tps_ord <- tryCatch(levels(order_tp(tps_all)), error = function(e) tps_all)
#   
#   # full grid to include zeros
#   regs <- c("upregulated", "downregulated")
#   all_grid <- tidyr::expand_grid(timepoint = tps_ord, regulation = regs)
#   
#   # counts + totals
#   counts <- df %>%
#     dplyr::count(timepoint, regulation, name = "n") %>%
#     dplyr::right_join(all_grid, by = c("timepoint", "regulation")) %>%
#     dplyr::mutate(n = tidyr::replace_na(n, 0L),
#                   timepoint = factor(timepoint, levels = tps_ord)) %>%
#     tidyr::pivot_wider(names_from = regulation, values_from = n, values_fill = 0) %>%
#     dplyr::arrange(timepoint) %>%
#     dplyr::mutate(Total = upregulated + downregulated) %>%
#     dplyr::rename(Upregulated = upregulated, Downregulated = downregulated) %>%
#     dplyr::mutate(timepoint = as.character(timepoint)) %>%
#     as.data.frame(stringsAsFactors = FALSE)
#   
#   if (return_data) return(counts)
#   
#   # Build caption
#   if (is.null(caption)) {
#     ent <- if ("gene" %in% names(TT)) "genes" else "entities"
#     caption <- sprintf(
#       "Dose-responsive %s per timepoint (%s < %s)",
#       ent, p_col, format(p_thresh)
#     )
#   }
#   
#   # Display table
#   tbl <- counts
#   names(tbl)[names(tbl) == "timepoint"] <- "Timepoint"
#   
#   tbl_output <- knitr::kable(
#     tbl,
#     caption = caption,
#     align = "l"
#   )
#   
#   invisible(counts)
#   
#   tbl_output
# 
# }




