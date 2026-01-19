plot_bmd_genes_top_per_tp <- function(
    top_pw_pertp,        # cols: timepoint, pathway_id, bmd, bmdl, bmdu  (pathway medians)
    gene_bmd,            # cols: timepoint, gene_id, bmd, bmdl, bmdu
    pathway_overview,    # cols: gene_id, pathway_id
    log10_x = TRUE,
    band_alpha = 0.07,
    band_fill = "pink",
    median_line_col = "red",
    show_in_knit = FALSE,
    label_threshold = 30,
    combine_per_tp = TRUE,
    ncol = 2,
    order_by = c("timepoint", "pathway_id")
) {
  stopifnot(
    is.data.frame(top_pw_pertp),
    is.data.frame(gene_bmd),
    is.data.frame(pathway_overview),
    all(c("timepoint","pathway_id","bmd","bmdl","bmdu") %in% names(top_pw_pertp)),
    all(c("timepoint","gene_id","bmd","bmdl","bmdu") %in% names(gene_bmd)),
    all(c("gene_id","pathway_id") %in% names(pathway_overview))
  )
  order_by <- match.arg(order_by, several.ok = TRUE)
  
  # Normalize keys (trim + character) to avoid silent mismatches
  norm <- function(x) trimws(as.character(x))
  top_pw_pertp <- top_pw_pertp |>
    dplyr::mutate(timepoint = norm(timepoint), pathway = norm(pathway_id))
  gene_bmd <- gene_bmd |>
    dplyr::mutate(timepoint = norm(timepoint), gene_id = norm(gene_id))
  pathway_overview <- pathway_overview |>
    dplyr::mutate(pathway = norm(pathway_id), gene_id = norm(gene_id))
  
  # ---- timepoint ordering helper (numeric) ----
  tp_num <- function(x) suppressWarnings(as.numeric(sub(".*?(\\d+).*", "\\1", x)))
  
  # Sort rows for stable output: numeric timepoint, then pathway
  pw_rows <- top_pw_pertp |>
    dplyr::distinct(timepoint, pathway, bmd, bmdl, bmdu) |>
    dplyr::mutate(.tp_num = tp_num(timepoint)) |>
    dplyr::arrange(.tp_num, timepoint, pathway) |>
    dplyr::select(-.tp_num)
  
  plots <- list()
  diag  <- dplyr::tibble(
    timepoint = character(), pathway = character(),
    n_genes_in_pathway = integer(), n_genes_with_bmd = integer(),
    note = character()
  )
  
  for (i in seq_len(nrow(pw_rows))) {
    tp  <- pw_rows$timepoint[i]
    pw  <- pw_rows$pathway[i]
    pwb <- pw_rows$bmd[i]; pwl <- pw_rows$bmdl[i]; pwu <- pw_rows$bmdu[i]
    
    genes_in_pw <- pathway_overview |>
      dplyr::filter(pathway == pw) |>
      dplyr::pull(gene_id) |>
      unique()
    
    if (length(genes_in_pw) == 0) {
      diag <- dplyr::add_row(diag, timepoint = tp, pathway = pw,
                             n_genes_in_pathway = 0, n_genes_with_bmd = 0,
                             note = "No genes mapped to this pathway in pathway_overview")
      next
    }
    
    dat <- gene_bmd |>
      dplyr::filter(timepoint == tp, gene_id %in% genes_in_pw) |>
      dplyr::mutate(gene_id = forcats::fct_reorder(gene_id, bmd))
    
    if (!nrow(dat)) {
      diag <- dplyr::add_row(diag, timepoint = tp, pathway = pw,
                             n_genes_in_pathway = length(genes_in_pw),
                             n_genes_with_bmd = 0,
                             note = "No BMD rows for these genes at this timepoint")
      next
    }
    
    p <- ggplot2::ggplot(dat, ggplot2::aes(x = bmd, y = gene_id)) +
      ggplot2::geom_pointrange(ggplot2::aes(xmin = bmdl, xmax = bmdu)) +
      ggplot2::labs(
        y = "Gene",
        x = "BMD (with BMDL/BMDU)",
        title = glue::glue("{pw} — {tp}")
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(panel.grid.major.y = ggplot2::element_line(color = "lightgray", linewidth = 0.5, linetype = 2)) +
      ggplot2::annotate("rect", xmin = pwl, xmax = pwu, ymin = -Inf, ymax = Inf,
                        fill = band_fill, alpha = band_alpha) +
      ggplot2::geom_vline(xintercept = pwb, linetype = 4, color = median_line_col, linewidth = 0.7)
    
    if (log10_x) p <- p + ggplot2::scale_x_log10()
    
    n_genes <- nrow(dat)
    if (n_genes > label_threshold) {
      p <- p +
        ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank()) +
        ggplot2::labs(y = NULL)
    }
    
    key <- paste(tp, pw, sep = " | ")
    plots[[key]] <- p
    
    diag <- dplyr::add_row(diag, timepoint = tp, pathway = pw,
                           n_genes_in_pathway = length(genes_in_pw),
                           n_genes_with_bmd = n_genes,
                           note = ifelse(n_genes > label_threshold, "OK (hidden y labels)", "OK"))
    
    if (isTRUE(show_in_knit)) print(p)
  }
  
  if (length(plots) == 0) {
    message("[plot_bmd_genes_top_per_tp] No plots produced.")
  }
  
  combined <- list()
  if (combine_per_tp && length(plots)) {
    keys <- names(plots)
    tps  <- sub(" \\| .*", "", keys)
    split_keys <- split(keys, tps)
    
    # ---- iterate timepoints in numeric order ----
    tp_levels <- unique(pw_rows$timepoint)  # already ordered
    for (tp in tp_levels) {
      if (!tp %in% names(split_keys)) next
      
      these_keys <- split_keys[[tp]]
      
      # keep same pathway order as pw_rows (already ordered)
      ordered_keys <- these_keys[
        order(match(these_keys, paste(pw_rows$timepoint, pw_rows$pathway, sep = " | ")))
      ]
      
      patch <- patchwork::wrap_plots(
        plots[ordered_keys],
        ncol = ncol,
        guides = "collect",
        byrow = TRUE
      ) + patchwork::plot_annotation(title = tp)
      
      combined[[tp]] <- patch
      if (isTRUE(show_in_knit)) print(patch)
    }
  }
  
  # return object
  if (combine_per_tp) {
    return(combined)  # named by timepoint
  } else {
    return(plots)     # named by "timepoint | pathway"
  }
}
