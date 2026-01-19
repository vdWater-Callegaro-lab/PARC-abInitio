

#' Accumulation plot of BMD across timepoints
#'
#' Build an accumulation plot (cumulative count vs. BMD) per timepoint, with
#' optional interactive output using Plotly that shows the gene/pathway label on hover.
#'
#' @param df A data.frame with at least columns `timepoint` and `bmd`. An
#'   optional label column (e.g. `gene`, `pathway`, `label`, `name`, `term`,
#'   `feature`) can be surfaced in hover tooltips.
#' @param title_prefix Character scalar used as the plot title prefix. Default "Gene".
#' @param log10_x Logical; if TRUE, the x-axis (BMD) is on a log10 scale.
#' @param palette Optional character vector of colors (hex or names) used for
#'   the timepoint curves. If NULL, a distinct qualitative palette is generated.
#' @param label_col Optional name of the column with a label to show on hover
#'   (e.g. a gene or pathway identifier). If NULL, the function will try common
#'   candidates: `c("gene","pathway","label","name","term","feature")`.
#' @param engine Output type: either "plotly" (default, interactive) or
#'   "ggplot" (static ggplot object).
#'
#' @return A Plotly object (htmlwidget) if `engine = "plotly"`, otherwise a ggplot object.
#'
#' @examples
#' # Basic usage (Plotly):
#' # p <- accumulation_plot(df)
#' # htmlwidgets::saveWidget(p, "accumulation.html")
#'
#' # Static ggplot:
#' # g <- accumulation_plot(df, engine = "ggplot")
#' # print(g)
#'
#' @importFrom dplyr group_by arrange mutate row_number ungroup filter
#' @importFrom ggplot2 ggplot aes geom_step geom_point scale_x_continuous
#' @importFrom ggplot2 scale_x_log10 scale_y_continuous scale_color_manual labs
#' @importFrom ggplot2 theme_minimal theme
#' @importFrom scales hue_pal log_breaks pretty_breaks
#' @importFrom glue glue
#' @importFrom stats setNames
#' @export
accumulation_plot <- function(
    df,
    title_prefix = "Gene",
    log10_x = FALSE,
    palette = NULL,
    label_col = NULL,
    engine = c("plotly", "ggplot")
) {
  engine <- match.arg(engine)
  stopifnot(is.data.frame(df), all(c("timepoint", "bmd") %in% names(df)))
  
  # -- 1) Clean / coerce -----------------------------------------------------------------
  d <- df %>%
    dplyr::mutate(bmd = as.numeric(.data$bmd)) %>%
    dplyr::filter(is.finite(.data$bmd))
  
  # -- 2) Natural ordering of timepoints (e.g., 4hr, 8hr, 1d) ---------------------------
  tps <- unique(as.character(d$timepoint))
  tps_ord <- tryCatch(.order_timepoints(tps), error = function(e) tps)
  d$timepoint <- factor(as.character(d$timepoint), levels = tps_ord)
  
  # -- 3) Cumulative ranks per timepoint -------------------------------------------------
  d_ranked <- d %>%
    dplyr::group_by(.data$timepoint) %>%
    dplyr::arrange(.data$bmd, .by_group = TRUE) %>%
    dplyr::mutate(rank = dplyr::row_number()) %>%
    dplyr::ungroup()
  
  # -- 4) Label column for hover ---------------------------------------------------------
  if (is.null(label_col)) {
    candidates <- c("gene","pathway_id","label","name","term","feature")
    label_col <- candidates[candidates %in% names(d_ranked)][1]
    if (length(label_col) == 0) label_col <- NULL
  } else {
    if (!label_col %in% names(d_ranked)) {
      warning(sprintf("label_col '%s' not found; hover will omit labels.", label_col))
      label_col <- NULL
    }
  }
  
  # -- 5) Colors ------------------------------------------------------------------------
  if (is.null(palette)) {
    n_colors <- length(levels(d_ranked$timepoint))
    palette <- scales::hue_pal()(n_colors)
  }
  # name colors to match timepoint order
  palette_named <- stats::setNames(palette, levels(d_ranked$timepoint))
  
  # -- 6) Build plot ---------------------------------------------------------------------
  title_text <- glue::glue("{title_prefix} BMD accumulation across timepoints")
  
  if (engine == "ggplot") {
    # ggplot (static)
    g <- ggplot2::ggplot(d_ranked, ggplot2::aes(x = .data$bmd, y = .data$rank, color = .data$timepoint)) +
      ggplot2::geom_step(linewidth = 0.25) +
      ggplot2::geom_point(size = 1, alpha = 0.6) +
      {
        if (log10_x) ggplot2::scale_x_log10("BMD", breaks = scales::log_breaks()) else ggplot2::scale_x_continuous("BMD", breaks = scales::pretty_breaks())
      } +
      ggplot2::scale_y_continuous("Cumulative count", breaks = scales::pretty_breaks()) +
      ggplot2::scale_color_manual(values = palette_named, name = "Timepoint") +
      ggplot2::labs(title = title_text, subtitle = "Each curve is a timepoint") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(legend.position = "right")
    return(g)
  }
  
  # Plotly (interactive)
  # Prepare hover text vector per-row
  if (!is.null(label_col)) {
    d_ranked$..label.. <- as.character(d_ranked[[label_col]])
  } else {
    d_ranked$..label.. <- NA_character_
  }
  
  plt <- plotly::plot_ly()
  
  tps_lvls <- levels(d_ranked$timepoint)
  for (tp in tps_lvls) {
    s <- d_ranked[d_ranked$timepoint == tp, , drop = FALSE]
    # Build hovertemplate: include label if available
    if (!all(is.na(s$..label..))) {
      hovertemplate <- paste0(
        "%{text}",
        "<br>Timepoint: ", tp,
        "<br>BMD: %{x}",
        "<br>Cumulative count: %{y}",
        "<extra></extra>"
      )
      text_vec <- s$..label..
    } else {
      hovertemplate <- paste0(
        "Timepoint: ", tp,
        "<br>BMD: %{x}",
        "<br>Cumulative count: %{y}",
        "<extra></extra>"
      )
      text_vec <- NULL
    }
    
    plt <- plt %>%
      plotly::add_trace(
        data = s,
        x = ~bmd, y = ~rank,
        type = "scatter", mode = "lines+markers",
        name = as.character(tp),
        line = list(shape = "hv", width = 1, color = palette_named[[as.character(tp)]]),
        marker = list(size = 6, opacity = 0.6, color = palette_named[[as.character(tp)]]),
        text = text_vec, hovertemplate = hovertemplate
      )
  }
  
  xaxis <- list(title = "BMD", hoverformat = ".2f")
  if (isTRUE(log10_x)) xaxis$type <- "log"
  
  plt <- plt %>%
    plotly::layout(
      title = list(text = title_text),
      xaxis = xaxis,
      yaxis = list(title = "Cumulative count"),
      legend = list(orientation = "v", x = 1.02, xanchor = "left", y = 1),
      hovermode = "closest"
    )
  
  plt
}

# ---- Helpers ---------------------------------------------------------------------------

#' @keywords internal
.order_timepoints <- function(x) {
  # Convert a character vector of timepoint labels into a sensible order.
  # Handles forms like: "4h", "4hr", "8hr", "24h", "1d", "2day", "3w", etc.
  stopifnot(is.character(x))
  if (!length(x)) return(x)
  
  # normalize
  z <- trimws(x)
  # extract numeric and unit
  num <- suppressWarnings(as.numeric(sub("^\\s*([0-9]*\\.?[0-9]+).*", "\\1", z)))
  unit_raw <- tolower(sub("^\\s*[0-9]*\\.?[0-9]+\\s*", "", z))
  
  # map units to minutes multiplier
  unit_map <- c(h = 60, hr = 60, hrs = 60, hour = 60, hours = 60,
                d = 1440, day = 1440, days = 1440,
                w = 10080, wk = 10080, week = 10080, weeks = 10080,
                min = 1, mins = 1, minute = 1, minutes = 1)
  
  # pick first recognized token from unit_raw
  unit_tok <- sub("^([a-zA-Z]+).*", "\\1", unit_raw)
  mult <- unit_map[unit_tok]
  mult[is.na(mult)] <- 60  # default to hours when unknown
  
  value_min <- num * as.numeric(mult)
  # If non-numeric, keep original order
  if (any(is.na(value_min))) return(x)
  
  ord <- order(value_min, na.last = TRUE)
  x[ord]
}


# accumulation_plot <- function(
#     df,
#     title_prefix = "Gene",
#     log10_x = FALSE,
#     palette = NULL
# ) {
#   stopifnot(is.data.frame(df), all(c("timepoint", "bmd") %in% names(df)))
#   
#   d <- df %>%
#     dplyr::mutate(bmd = as.numeric(bmd)) %>%
#     dplyr::filter(is.finite(bmd))
#   
#   # Natural ordering of timepoints (e.g., 4hr, 8hr, 24hr)
#   tps <- unique(as.character(d$timepoint))
#   tps_ord <- tryCatch(levels(order_tp(tps)), error = function(e) tps)
#   d$timepoint <- factor(as.character(d$timepoint), levels = tps_ord)
#   
#   # Cumulative ranks per timepoint
#   d_ranked <- d %>%
#     dplyr::group_by(timepoint) %>%
#     dplyr::arrange(bmd, .by_group = TRUE) %>%
#     dplyr::mutate(rank = dplyr::row_number()) %>%
#     dplyr::ungroup()
#   
#   # Colors
#   if (is.null(palette)) {
#     n_colors <- length(unique(d_ranked$timepoint))
#     palette <- scales::hue_pal()(n_colors)
#   }
#   
#   ggplot(d_ranked, aes(x = bmd, y = rank, color = timepoint)) +
#     geom_step(linewidth = 0.2) +
#     geom_point(size = 0.8, alpha = 0.6) +
#     { if (log10_x) scale_x_log10("BMD", breaks = scales::log_breaks()) else scale_x_continuous("BMD", breaks = scales::pretty_breaks()) } +
#     scale_y_continuous("Cumulative count", breaks = scales::pretty_breaks()) +
#     scale_color_manual(values = palette, name = "Timepoint") +
#     labs(
#       title = glue::glue("{title_prefix} BMD accumulation across timepoints"),
#       subtitle = "Each curve is a timepoint"
#     ) +
#     theme_minimal(base_size = 11) +
#     theme(legend.position = "right")
# }
