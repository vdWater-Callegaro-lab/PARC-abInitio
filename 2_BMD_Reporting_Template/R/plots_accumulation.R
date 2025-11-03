

accumulation_plot <- function(
    df,
    title_prefix = "Gene",
    log10_x = FALSE,
    palette = NULL
) {
  stopifnot(is.data.frame(df), all(c("timepoint", "bmd") %in% names(df)))
  
  d <- df %>%
    dplyr::mutate(bmd = as.numeric(bmd)) %>%
    dplyr::filter(is.finite(bmd))
  
  # Natural ordering of timepoints (e.g., 4hr, 8hr, 24hr)
  tps <- unique(as.character(d$timepoint))
  tps_ord <- tryCatch(levels(order_tp(tps)), error = function(e) tps)
  d$timepoint <- factor(as.character(d$timepoint), levels = tps_ord)
  
  # Cumulative ranks per timepoint
  d_ranked <- d %>%
    dplyr::group_by(timepoint) %>%
    dplyr::arrange(bmd, .by_group = TRUE) %>%
    dplyr::mutate(rank = dplyr::row_number()) %>%
    dplyr::ungroup()
  
  # Colors
  if (is.null(palette)) {
    n_colors <- length(unique(d_ranked$timepoint))
    palette <- scales::hue_pal()(n_colors)
  }
  
  ggplot(d_ranked, aes(x = bmd, y = rank, color = timepoint)) +
    geom_step(linewidth = 1) +
    geom_point(size = 0.8, alpha = 0.6) +
    { if (log10_x) scale_x_log10("BMD", breaks = scales::log_breaks()) else scale_x_continuous("BMD", breaks = scales::pretty_breaks()) } +
    scale_y_continuous("Cumulative count", breaks = scales::pretty_breaks()) +
    scale_color_manual(values = palette, name = "Timepoint") +
    labs(
      title = glue::glue("{title_prefix} BMD accumulation across timepoints"),
      subtitle = "Each curve is a timepoint"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "right")
}
