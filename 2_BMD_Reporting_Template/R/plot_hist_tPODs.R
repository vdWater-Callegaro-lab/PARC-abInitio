
# Histogram of BMD values with tPOD reference lines
plot_bmd_hist_with_tpods <- function(
    bmd,                               # numeric vector OR data.frame with columns value_col and optional tp_col
    tpods,                             # data.frame with tPOD metrics; optional tp_col to match facets
    value_col   = "bmd",               # column name in bmd (ignored if bmd is numeric)
    tp_col      = NULL,                # timepoint column name in bmd (for faceting). e.g., "timepoint"
    tpods_tpcol = "timepoint",         # timepoint column name in tpods (if present)
    # name mapping of your tPOD columns:
    metric_cols = c(P5 = "P5", Nth25 = "Nth25", LCRD = "LCRD", FirstMode = "FirstMode", Kneedle = "Kneedle"),
    bins        = 30,
    density     = FALSE,               # add density outline on top of histogram
    log10_x     = FALSE,               # log10 axis for BMDs
    palette     = c(P5 = "#1f77b4", Nth25 = "#2ca02c", LCRD = "#d62728", FirstMode = "#9467bd", Kneedle = "#ff7f0e"),
    title       = NULL,
    subtitle    = NULL,
    band_alpha  = 0.0,                 # set >0 to add a faint band around each tPOD? (kept off by default)
    show_in_knit = FALSE               # print immediately (useful in simple chunks; inside if/loops still prefer print(p))
) {
  stopifnot(is.data.frame(tpods) || is.list(tpods))
  
  # Coerce BMD input to a data.frame with columns: value_col and optional tp_col
  if (is.numeric(bmd)) {
    bmd_df <- data.frame(!!value_col := as.numeric(bmd))
  } else if (is.data.frame(bmd)) {
    if (!value_col %in% names(bmd)) stop(sprintf("Column '%s' not found in 'bmd'.", value_col))
    bmd_df <- bmd
  } else stop("`bmd` must be a numeric vector or a data.frame.")
  
  # Clean BMD values
  bmd_df <- bmd_df[is.finite(bmd_df[[value_col]]), , drop = FALSE]
  
  # Verify tPOD metric columns exist (collect only those that are present)
  metric_cols_present <- metric_cols[metric_cols %in% names(tpods)]
  if (length(metric_cols_present) == 0L) {
    stop(
      "None of the specified metric columns were found in `tpods`. ",
      "Looked for: ", paste(unname(metric_cols), collapse = ", "),
      " — found in tpods: ", paste(names(tpods), collapse = ", ")
    )
  }
  
  # Long-format tPODs
  # IMPORTANT: use UNNAMED column vector for `cols=` to avoid the rename error
  cols_for_pivot <- unname(metric_cols_present)
  
  tpods_long <- tidyr::pivot_longer(
    tpods,
    cols = tidyselect::all_of(cols_for_pivot),
    names_to = "metric_col",
    values_to = "value"
  )
  
  # Map column names -> friendly metric labels (p5, rank25, lcrd, first_mode)
  # metric_cols_present is still NAMED (names = friendly labels, values = column names)
  colname_to_label <- stats::setNames(names(metric_cols_present), unname(metric_cols_present))
  tpods_long$metric <- colname_to_label[tpods_long$metric_col]
  
  # keep only finite values
  tpods_long <- tpods_long[is.finite(tpods_long$value), , drop = FALSE]
  
  # If faceting by timepoint: ensure both have the same tp col; else treat as global (no facet)
  facetting <- !is.null(tp_col) && tp_col %in% names(bmd_df)
  if (facetting) {
    # keep only tpods columns needed
    if (!(tpods_tpcol %in% names(tpods_long))) {
      stop(sprintf("Faceting requested by '%s', but '%s' not found in `tpods`.", tp_col, tpods_tpcol))
    }
    # coerce to character to avoid level mismatches
    bmd_df[[tp_col]]     <- as.character(bmd_df[[tp_col]])
    tpods_long[[tpods_tpcol]] <- as.character(tpods_long[[tpods_tpcol]])
    # filter tpods to timepoints present in bmd_df
    tpods_long <- tpods_long[tpods_long[[tpods_tpcol]] %in% unique(bmd_df[[tp_col]]), , drop = FALSE]
  }
  
  # Base plot
  p <- ggplot(bmd_df, aes(x = .data[[value_col]])) +
    geom_histogram(bins = bins, fill = "grey80", color = "grey40")
  
  if (density) {
    p <- p + 
      ggplot2::geom_density(bw = "nrd0", linewidth = 0.7, alpha = 0.0)
  }
  
  # Add tPOD lines (with legend)
  if (facetting) {
    vline_df <- tpods_long |>
      dplyr::rename(!!tp_col := !!tpods_tpcol)
  } else {
    vline_df <- tpods_long
  }
  
  if (nrow(vline_df)) {
    p <- p +
      geom_vline(
        data = vline_df,
        aes(xintercept = value, color = metric),
        linetype = 2, linewidth = 0.9, show.legend = TRUE
      )
    # optional faint band around the line (off by default)
    if (band_alpha > 0) {
      p <- p +
        annotate("rect",
                 xmin = vline_df$value * 0.999, xmax = vline_df$value * 1.001,
                 ymin = -Inf, ymax = Inf,
                 alpha = band_alpha)
    }
  }
  
  # Facet by timepoint if asked
  if (facetting) {
    p <- p + facet_wrap(stats::as.formula(paste("~", tp_col)), scales = "free_y")
  }
  
  # Scales & labels
  p <- p +
    scale_color_manual(
      name = "tPOD",
      breaks = names(palette)[names(palette) %in% unique(vline_df$metric)],
      values = palette
    ) +
    labs(
      x = "BMD",
      y = "Count",
      title = title %||% "BMD histogram with tPOD overlays",
      subtitle = subtitle %||% if (facetting) NULL else NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  if (log10_x) {
    p <- p + scale_x_log10()
  }
  
  if (isTRUE(show_in_knit)) print(p)
  return(p)
}

# Helper: safe `%||%`
`%||%` <- function(a, b) if (!is.null(a)) a else b
