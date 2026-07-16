
# functions


# BMD accumulation plot
BMDaccumulationPlot <- function(BMDoutput, analysis_summary) {
  BMDoutput %>%
    group_by(timepoint) %>%
    arrange(finalBMD) %>%
    mutate(ranked_gene = row_number()) %>%
    ggplot(aes(finalBMD, y = ranked_gene)) +
    geom_point(stat = "identity", aes(color = timepoint, shape = timepoint), size = 0.8) +
    labs(title = analysis_summary,
         y = "Accumulation",
         x = "BMC") +
    # scale_color_manual(values = c("#E6E6E6", "#D1D1D1","#BBBBBB","#A0A0A0","#7F7F7F","#4D4D4D"))
    scale_color_manual(values = timepoint_cols)
}




# plot BMD accumulation genes zoomed in
BMDaccumulationPlot_zoom <- function(BMDoutput) {
  BMDoutput %>%
    group_by(timepoint) %>%
    arrange(finalBMD) %>%
    mutate(ranked_gene = row_number()) %>%
    ggplot(aes(finalBMD, y = ranked_gene)) +
    geom_point(stat = "identity", aes(color = timepoint, shape = timepoint), size = 1) +
    labs(y = "Accumulation", 
         x = "BMC") +
    ylim(c(0, 250)) + 
    scale_color_manual(values = timepoint_cols)
}



# plot median BMD accumulation of pathways
BMDaccumulationPlot_pathway <- function(BMDoutput, analysis_summary) {
  BMDoutput %>%
    group_by(timepoint) %>%
    arrange(medianBMD) %>%
    mutate(ranked_pathway = row_number()) %>%
    ggplot(aes(medianBMD, y = ranked_pathway)) +
    geom_point(stat = "identity", aes(color = timepoint, shape = timepoint), size = 1.5) +
    labs(title = analysis_summary,
         y = "Accumulation", 
         x = "Median BMC") +
    # scale_color_manual(values = c("#E6E6E6", "#D1D1D1","#BBBBBB","#A0A0A0","#7F7F7F","#4D4D4D"),
    #                    breaks = c("4h", "8h", "16h", "24h", "48h", "72h")) +
    scale_color_manual(values = c(timepoint_cols),
                       breaks = c("4h", "8h", "16h", "24h", "48h", "72h")) +
    scale_shape_manual(values = c(16, 17, 15, 3, 7, 8),
                       breaks = c("4h", "8h", "16h", "24h", "48h", "72h"))
}




# plot median BMD accumulation of pathways + highlight pathway_oi
BMDaccumulationPlot_pathway_highlight <- function(BMDoutput, analysis_summary, pathway_oi) {
  plot_data =  BMDoutput %>%
    group_by(timepoint) %>%
    arrange(medianBMD) %>%
    mutate(ranked_pathway = row_number())

  highlight_data = plot_data %>% filter(Pathway.Name %in% pathway_oi)

    ggplot(plot_data, aes(x= medianBMD, y = ranked_pathway)) +
    geom_point(stat = "identity", aes(color = timepoint, shape = timepoint), size = 1.5) +
    geom_point(data = highlight_data, aes(x = medianBMD, y = ranked_pathway, shape = timepoint), color = "red", size = 3) +
    # geom_text(data = highlight_data, aes(x = medianBMD, y = ranked_pathway, label = Pathway.Name), vjust = -1, hjust = 1, color = "red") +
    labs(title = analysis_summary,
         y = "Accumulation",
         x = "Median BMC") +
    # scale_color_manual(values = c("#E6E6E6", "#D1D1D1","#BBBBBB","#A0A0A0","#7F7F7F","#4D4D4D"))
      scale_color_manual(values = timepoint_cols)

}




# plot median BMD accumulation of pathways for all partners combined in one plot + highlight pathway_oi
BMDaccumulationPlot_pathway_highlight_combined <- function(BMDoutput, timepoint_oi, pathway_oi) {
  plot_data = BMDoutput %>%
    filter(timepoint == timepoint_oi) %>%
    group_by(analysis_summary) %>%
    arrange(medianBMD) %>%
    mutate(ranked_pathway = row_number()) 
  
  highlight_data = plot_data %>% filter(Pathway.Name %in% pathway_oi)
  
  ggplot(plot_data, aes(x= medianBMD, y = ranked_pathway)) +
    geom_point(stat = "identity", aes(color = analysis_summary), size = 1.5) +
    geom_point(data = highlight_data, aes(x = medianBMD, y = ranked_pathway, color = analysis_summary), size = 5) +
    labs(y = "Accumulation",
         x = "Median BMC",
         subtitle = timepoint_oi) +
    scale_color_manual(values = c("blue2", "forestgreen","firebrick", "gray30", "gold1"),
                       breaks = c( "BMDE-noWTT-CPM-RF-S5", "BMDE-WTT-CPM-RF-S0", "DRO-Quad-UQ-RF-S0", "DRO-Quad-VST-C10-S0", "DRO-Quad-VST-RF-S0")) +
    scale_shape_manual(values = c(15, 16, 17, 18, 8)) +
    theme(legend.title = element_blank())
  
}





# plot BMD of genes in a pathway
plot_BMD_genes_pathwayoi = function(genes_in_pathway, medianBMD_pathways, time_oi, pathway_oi, BMDoutput_genes, pathway_overview) {
    
    BMDL_pathway = medianBMD_pathways %>% filter(timepoint == time_oi & Pathway.Name == pathway_oi) %>% pull(medianBMDL) %>% unique()
    BMD_pathway = medianBMD_pathways %>% filter(timepoint == time_oi & Pathway.Name == pathway_oi) %>% pull(medianBMD) %>% unique()
    BMDU_pathway = medianBMD_pathways %>% filter(timepoint == time_oi & Pathway.Name == pathway_oi) %>% pull(medianBMDU) %>% unique()
    
    plot_data = BMDoutput_genes %>% separate(gene_symbol_entrez_id, into = c("gene_symbol", "entrez_id"), sep = "_") %>% 
      filter(timepoint == time_oi & gene_symbol %in% genes_in_pathway)
    
    p = ggplot(plot_data, aes(x = finalBMD, y = reorder(gene_symbol, finalBMD))
    ) +
      geom_pointrange(aes(xmin = finalBMDL, xmax = finalBMDU), fatten = 1) +
      scale_x_log10() +
      labs(y = "Gene Symbol",
           x = "BMD (BMDL/BMDU)") +
      geom_rect(
        aes(
          xmin = BMDL_pathway,
          xmax = BMDU_pathway,
          ymin = -Inf,
          ymax = Inf
        ),
        fill = "pink",
        alpha = 0.01
      ) +
      geom_vline(
        xintercept = BMD_pathway,
        lty = 4,
        color = "red",
        linewidth = 0.7
      ) +
      theme(panel.grid.major.y = element_line(
        color = "lightgray",
        linewidth = 0.5,
        linetype = 2
      ),
      axis.text.y = element_text(size = 6),
      axis.text.x = element_text(size = 8))
    
    return(p)
}




# BMD accumulation plot for BMD ran on pw scores
BMDaccumulationPlot_pwscores_highlight <- function(BMDoutput, HLpw, analysis_method){
  plotdata = BMDoutput %>%
    group_by(timepoint) %>%
    arrange(finalBMD) %>%
    mutate(ranked_pathway = row_number())
  
  highlight_data = plotdata %>% filter(pathway == HLpw)
  
  ggplot(plotdata, aes(finalBMD, y = ranked_pathway)) + 
    geom_point(stat = "identity", aes(color = timepoint, shape = timepoint), size = 1.5) +
    geom_point(data = highlight_data, aes(x = finalBMD, y=ranked_pathway, shape = timepoint), color = "red", size = 3) +
    labs(title = analysis_method,
         y = "Accumulation", 
         x = "BMC") +
    scale_color_manual(values = c(timepoint_cols))

}




# tPOD CALCULATIONS
# 5th percentile (or any percentile)
get_percentile_bmd <- function(df, prob = 0.05, na.rm = TRUE) {
  df %>%
    group_by(timepoint) %>%
    summarise(tpod = quantile(bmd, probs = prob, na.rm = na.rm),
              .groups = "drop")
}

# N-th ranked (e.g., 25th)
get_nth_ranked_bmd <- function(df, n = 25) {
  df %>%
    group_by(timepoint) %>%
    arrange(bmd, .by_group = TRUE) %>%
    mutate(row = row_number()) %>%
    filter(row == n) %>%
    summarise(tpod = first(bmd), .groups = "drop")
}

# Lowest Consistent Ranked Dose (LCRD)
get_LCRD_bmd <- function(df, ratio_threshold = 1.66) {
  compute_lcrd <- function(x) {
    x <- sort(x)
    if (length(x) < 2) return(NA_real_)
    groups <- list(current = x[1])
    largest <- groups$current
    for (i in 2:length(x)) {
      if (x[i] / x[i - 1] <= ratio_threshold) {
        groups$current <- c(groups$current, x[i])
      } else {
        if (length(groups$current) > length(largest)) largest <- groups$current
        groups$current <- x[i]
      }
    }
    if (length(groups$current) > length(largest)) largest <- groups$current
    min(largest)
  }
  
  df %>%
    group_by(timepoint) %>%
    summarise(tpod = compute_lcrd(bmd), .groups = "drop")
}

# First mode of the BMD distribution
get_first_mode_bmd <- function(df, bw = "nrd0") {
  first_mode <- function(x) {
    x <- x[is.finite(x)]
    if (length(unique(x)) < 2) return(NA_real_)
    dens <- density(x, bw = bw)
    y <- dens$y
    xgrid <- dens$x
    idx <- which(diff(sign(diff(y))) == -2) + 1  # local maxima
    if (length(idx) == 0) return(NA_real_)
    xgrid[min(idx)]  # first mode
  }
  
  df %>%
    group_by(timepoint) %>%
    summarise(tpod = first_mode(bmd), .groups = "drop")
}

# Kneedle / inflection-point BMD (per timepoint), matching LCRD structure
get_kneedle_bmd <- function(df, S = 1, min_points = 3) {
  stopifnot(is.data.frame(df), all(c("timepoint", "bmd") %in% names(df)))
  
  compute_knee <- function(y_raw, S = 1, min_points = 3) {
    # need enough points to detect an inflection
    y_raw <- y_raw[is.finite(y_raw)]
    if (length(y_raw) < min_points) return(NA_real_)
    
    # Order by BMD; x is rank, y is sorted BMD
    ord <- order(y_raw)
    x   <- seq_along(y_raw)
    y   <- y_raw[ord]
    
    # Smooth the curve (fall back to raw if spline fails)
    ss <- tryCatch(stats::smooth.spline(x, y), error = function(e) NULL)
    if (!is.null(ss)) { x_s <- ss$x; y_s <- ss$y } else { x_s <- x; y_s <- y }
    
    # Normalize x,y to [0,1] (guard zero ranges)
    rx <- range(x_s, finite = TRUE); ry <- range(y_s, finite = TRUE)
    x_n <- if (diff(rx) == 0) rep(0, length(x_s)) else (x_s - rx[1]) / diff(rx)
    y_n <- if (diff(ry) == 0) rep(0, length(y_s)) else (y_s - ry[1]) / diff(ry)
    
    # Difference curve D(x) = y - x
    Dd <- data.frame(x = x_n, y = y_n - x_n, i = seq_along(x_n))
    
    # Find local maxima of D(x)
    if (nrow(Dd) < 3) {
      knee_i <- which.max(Dd$y)
    } else {
      lmx_idx <- which(diff(sign(diff(Dd$y))) == -2) + 1
      if (length(lmx_idx) == 0) {
        knee_i <- which.max(Dd$y)
      } else {
        Dlmx   <- Dd[lmx_idx, , drop = FALSE]
        diff_x <- mean(diff(Dd$x))
        Tlmx   <- Dlmx$y - S * diff_x
        
        # Search regions after each local max for first index below its threshold
        knee_idx <- integer(0)
        n <- nrow(Dd)
        for (k in seq_len(nrow(Dlmx))) {
          start <- Dlmx$i[k]
          end   <- if (k < nrow(Dlmx)) Dlmx$i[k + 1] else n
          idx   <- start:end
          below <- which(Dd$y[idx] < Tlmx[k])
          if (length(below) > 0) knee_idx <- c(knee_idx, idx[below])
        }
        
        knee_i <- if (length(knee_idx) == 0) which.max(Dd$y) else min(knee_idx)
      }
    }
    
    # Map knee index back to the ORIGINAL (unsmoothed) sorted y
    as.numeric(y[knee_i])
  }
  
  df %>%
    dplyr::group_by(timepoint) %>%
    dplyr::summarise(tpod = compute_knee(bmd, S = S, min_points = min_points),
                     .groups = "drop")
}




# Summarize BMD metrics per timepoint (simple text table; no knitr/kableExtra)
bmd_summary_table <- function(
    df,
    digits = 3,
    caption = "BMD summary per timepoint",
    return_data = FALSE
) {
  stopifnot(is.data.frame(df), all(c("timepoint","bmd") %in% names(df)))
  
  # Base timepoint order (natural if order_tp is available)
  tps <- unique(as.character(df$timepoint))
  tps_ord <- tryCatch(levels(order_tp(tps)), error = function(e) tps)
  base <- tibble::tibble(timepoint = factor(tps_ord, levels = tps_ord))
  
  # Compute metrics (using your helpers)
  p5    <- get_percentile_bmd(df, prob = 0.05)        |> dplyr::rename(`P5` = tpod)
  n25   <- get_nth_ranked_bmd(df, n = 25)             |> dplyr::rename(`Nth25` = tpod)
  lcrd  <- get_LCRD_bmd(df)                           |> dplyr::rename(`LCRD` = tpod)
  mode1 <- get_first_mode_bmd(df)                     |> dplyr::rename(`FirstMode` = tpod)
  kneedle <- get_kneedle_bmd(df)                      |> dplyr::rename(`Kneedle` = tpod)
  
  # Join everything together (keeping all timepoints)
  out <- base |>
    dplyr::left_join(p5,    by = "timepoint") |>
    dplyr::left_join(n25,   by = "timepoint") |>
    dplyr::left_join(lcrd,  by = "timepoint") |>
    dplyr::left_join(mode1, by = "timepoint") |>
    dplyr::left_join(kneedle, by = "timepoint") |>
    dplyr::mutate(timepoint = as.character(timepoint))
  
  # Round numerics for display
  num_cols <- names(out)[sapply(out, is.numeric)]
  out[num_cols] <- lapply(out[num_cols], function(x) round(x, digits))
  
  if (return_data) return(out)
  
  # ---------- simple monospaced table printing ----------
  tbl <- out
  names(tbl)[names(tbl) == "timepoint"] <- "Timepoint"
  
  # Format numeric cols to fixed decimals; keep character cols as-is
  fmt_num <- function(x) ifelse(is.na(x), "NA", sprintf(paste0("%.", digits, "f"), x))
  for (nm in names(tbl)) {
    if (is.numeric(tbl[[nm]])) {
      tbl[[nm]] <- fmt_num(tbl[[nm]])
    } else {
      tbl[[nm]] <- as.character(tbl[[nm]])
    }
  }
  
  # Column widths: max of header vs values
  widths <- sapply(names(tbl), function(nm) {
    max(nchar(nm), if (nrow(tbl) > 0) max(nchar(tbl[[nm]]), na.rm = TRUE) else 0)
  })
  
  # Alignment: left for first column, right for the rest
  fmt_cell <- function(val, w, right = FALSE) {
    if (right) sprintf(paste0("%", w, "s"), val) else sprintf(paste0("%-", w, "s"), val)
  }
  align_right <- c(FALSE, rep(TRUE, length(widths) - 1))
  
  # Header + rows
  header <- paste(mapply(function(nm, w, r) fmt_cell(nm, w, r), names(tbl), widths, align_right), collapse = "  ")
  sep_line <- paste(vapply(widths, function(w) paste(rep("-", w), collapse = ""), character(1L)), collapse = "  ")
  rows <- apply(tbl, 1, function(r) {
    paste(mapply(function(val, w, r) fmt_cell(val, w, r), r, widths, align_right), collapse = "  ")
  })
  
  # Caption and table
  if (!is.null(caption) && nzchar(caption)) cat(caption, "\n")
  cat(header, "\n", sep_line, "\n", sep = "")
  if (length(rows)) cat(paste(rows, collapse = "\n"), "\n", sep = "")
  invisible(out)
}





# # plot dose response curves modules > now does not work anymore; can probably be removed as we do not have this in our figure anymore
# BMD_output_plots_pathways = function(BMDoutput, time_oi, BMD_pathways_oi, txgmapr_egs) {
#   
#   if((BMDoutput %>% filter(timepoint == time_oi) %>% nrow()) == 0){
#     print("There are no pathways showing a dose-response")
#   }else if((BMDoutput %>% filter(timepoint == time_oi) %>% nrow()) < 5){
#     print("NOTE: There are less than 5 pathways showing a dose-response")
#   } 
#   
#   for (pathway in BMD_pathways_oi) {
#     BMD_pathway = pathway
#     
#     BMD_data <-
#       BMDoutput %>% filter(patwhay == BMD_pathway &
#                                      timepoint == time_oi)
#     
#     txg_egs_meta_data = txgmapr_egs %>% as.data.frame() %>% # select module of interest from txg data
#       filter(module == BMD_module) %>%
#       separate(
#         sample_id,
#         into = c(
#           "upload",
#           "exposure_type",
#           "replicate", 
#           "timepoint", 
#           "dose"
#         ),
#         sep="_",
#         remove=FALSE
#       ) %>%
#       mutate(
#       #   dose = case_when(
#       #   dose_level == "1" ~ 0.1,
#       #   dose_level == "2" ~ 0.5,
#       #   dose_level == "3" ~ 1,
#       #   dose_level == "4" ~ 2.5,
#       #   dose_level == "5" ~ 5,
#       #   dose_level == "6" ~ 10,
#       #   dose_level == "7" ~ 20,
#       #   dose_level == "8" ~ 30,
#       #   dose_level == "9" ~ 50
#       # ),
#       timepoint = paste0(timepoint, "h"),
#       dose = as.numeric(dose),
#       replicate = gsub("CISPLATIN", "", replicate),
#       replicate = paste0("R", replicate)) %>%
#       filter(timepoint == time_oi) 
#     
#     
#     BMD_padj = BMD_data$Prefilter.Adjusted.P.Value
#     BMD_padj_gene = BMD_data$Min.p.adjust
#     # BMD_corEG = BMD_data$corEG
#     BMDL = BMD_data$finalBMDL
#     BMDU = BMD_data$finalBMDU
#     # Cmax_low = Cmax_data$Cmax_low[Cmax_data$COMPOUND == BMD_compound]
#     # Cmax_high = Cmax_data$Cmax_high[Cmax_data$COMPOUND == BMD_compound]
#     BMD_y = BMD_data$finalBMD
#     BMD_model = BMD_data$Best.Model
#     BMD_module = BMD_data$module
#     
#     # if (BMD_data$Max.Fold.Change > 0) {
#     #   BMD_logFC = log2(BMD_data$Max.Fold.Change.Absolute.Value)
#     # } else {
#     #   BMD_logFC = log2(1 / BMD_data$Max.Fold.Change.Absolute.Value)
#     # }
#     
#     HILL_INTERCEPT = BMD_data$Hill.Parameter.Intercept
#     HILL_V = BMD_data$Hill.Parameter.v
#     HILL_N = BMD_data$Hill.Parameter.n
#     HILL_K = BMD_data$Hill.Parameter.k
#     
#     POWER_CONTROL = BMD_data$Power.Parameter.control
#     POWER_SLOPE = BMD_data$Power.Parameter.slope
#     POWER_POWER = BMD_data$Power.Parameter.power
#     
#     LINEAR_0 = BMD_data$Linear.Parameter.beta_0
#     LINEAR_1 = BMD_data$Linear.Parameter.beta_1
#     
#     POLY2_0 = BMD_data$Poly.2.Parameter.beta_0
#     POLY2_1 = BMD_data$Poly.2.Parameter.beta_1
#     POLY2_2 = BMD_data$Poly.2.Parameter.beta_2
#     
#     EXP2_A = BMD_data$Exp.2.Parameter.a
#     EXP2_B = BMD_data$Exp.2.Parameter.b
#     EXP2_SIGN = BMD_data$Exp.2.Parameter.sign
#     
#     EXP3_A = BMD_data$Exp.3.Parameter.a
#     EXP3_B = BMD_data$Exp.3.Parameter.b
#     EXP3_D = BMD_data$Exp.3.Parameter.d
#     EXP3_SIGN = BMD_data$Exp.3.Parameter.sign
#     
#     EXP4_A = BMD_data$Exp.4.Parameter.a
#     EXP4_B = BMD_data$Exp.4.Parameter.b
#     EXP4_C = BMD_data$Exp.4.Parameter.c
#     
#     EXP5_A = BMD_data$Exp.5.Parameter.a
#     EXP5_B = BMD_data$Exp.5.Parameter.b
#     EXP5_C = BMD_data$Exp.5.Parameter.c
#     EXP5_D = BMD_data$Exp.5.Parameter.d
#     
#     
#     
#     HILL   = function(DOSE)
#       HILL_INTERCEPT + HILL_V * DOSE ^ HILL_N / (HILL_K ^ HILL_N + DOSE ^ HILL_N)
#     POWER  = function(DOSE)
#       POWER_CONTROL + POWER_SLOPE * DOSE ^ POWER_POWER
#     LINEAR = function(DOSE)
#       LINEAR_0 + LINEAR_1 * DOSE
#     POLY2  = function(DOSE)
#       POLY2_0 + POLY2_1 * DOSE + POLY2_2 * DOSE ^ 2
#     EXP2   = function(DOSE)
#       EXP2_A * exp(EXP2_SIGN * EXP2_B * DOSE)
#     EXP3   = function(DOSE)
#       EXP3_A * exp(EXP3_SIGN * (EXP3_B * DOSE) ^ EXP3_D)
#     EXP4   = function(DOSE)
#       EXP4_A * (EXP4_C - (EXP4_C - 1) * exp(-EXP4_B * DOSE))
#     EXP5   = function(DOSE)
#       EXP5_A * (EXP5_C - (EXP5_C - 1) * exp(-(EXP5_B * DOSE) ^ EXP5_D))
#     
#     if (BMD_model == "Hill") {
#       BMD_function = HILL
#     } else if (BMD_model == "Power") {
#       BMD_function = POWER
#     } else if (BMD_model == "Linear") {
#       BMD_function = LINEAR
#     } else if (BMD_model == "Poly 2") {
#       BMD_function = POLY2
#     } else if (BMD_model == "Exp 2") {
#       BMD_function = EXP2
#     } else if (BMD_model == "Exp 3") {
#       BMD_function = EXP3
#     } else if (BMD_model == "Exp 4") {
#       BMD_function = EXP4
#     } else if (BMD_model == "Exp 5") {
#       BMD_function = EXP5
#     }
#     
#     yaxis_lim = c(min(txg_egs_meta_data$eg_score),
#                   max(txg_egs_meta_data$eg_score))
#     
#     xaxis_lim = c(min(txg_egs_meta_data$dose),
#                   max(txg_egs_meta_data$dose))
#     
#     
#     p = ggplot(data = txg_egs_meta_data, mapping = aes(x = dose, y = eg_score)) +
#       # theme_standard() +
#       ggtitle(
#         paste0(time_oi, " - ", BMD_module, " - ", BMD_model),
#         subtitle = paste0("BMD: ", round(BMD_y, 2), " µM, WTT p-adj: ", round(BMD_padj, 5))
#       ) +
#       geom_point(aes(group = replicate, shape = replicate), size = 3) +
#       stat_function(fun = BMD_function,
#                     colour = "blue",
#                     size = 1) +
#       geom_rect(
#         aes(
#           xmin = BMDL,
#           xmax = BMDU,
#           ymin = -Inf,
#           ymax = Inf
#         ),
#         fill = "pink",
#         alpha = 0.01
#       ) +
#       ylim(yaxis_lim) +
#       xlab("Concentration (µM)") + ylab("Module eigengene score") +
#       scale_shape_manual(name = "Replicate",
#                          labels = c("R1", "R2", "R3"),
#                          values = rep(seq(21, 25), 5)) +
#       geom_vline(xintercept = BMD_y, lty=4, color="blue", linewidth=1) +
#       scale_x_log10() 
#     
#     return(p)
#     
#   }
# 
# }



# plot_BMD_genes_moduleBMD = function(BMDmodules, time_oi, module_oi, BMDoutput_genes, pathway_overview) {
#   
#   BMDL_pathway = BMDmodules %>% filter(timepoint == time_oi & module == module_oi) %>% pull(finalBMDL) %>% unique()
#   BMD_pathway = BMDmodules %>% filter(timepoint == time_oi & module == module_oi) %>% pull(finalBMD) %>% unique()
#   BMDU_pathway = BMDmodules %>% filter(timepoint == time_oi & module == module_oi) %>% pull(finalBMDU) %>% unique()
#   
#   genes_in_pathway = pathway_overview %>% filter(module == module_oi) %>% 
#     separate(gene_symbol_entrez_id, into = c("gene_symbol", "entrez_id"), sep = "_") %>% 
#     pull(gene_symbol) %>% unique()
#   
#   plot_data = BMDoutput_genes %>% separate(gene_symbol_entrez_id, into = c("gene_symbol", "entrez_id"), sep = "_") %>% 
#     filter(timepoint == time_oi & gene_symbol %in% genes_in_pathway)
#   
#   p = ggplot(plot_data, aes(x = finalBMD, y = reorder(gene_symbol, finalBMD))
#   ) +
#     geom_pointrange(aes(xmin = finalBMDL, xmax = finalBMDU), fatten = 1) +
#     scale_x_log10() +
#     labs(y = "Gene Symbol",
#          x = "BMD (BMDL/BMDU)") +
#     geom_rect(
#       aes(
#         xmin = BMDL_pathway,
#         xmax = BMDU_pathway,
#         ymin = -Inf,
#         ymax = Inf
#       ),
#       fill = "pink",
#       alpha = 0.01
#     ) +
#     geom_vline(
#       xintercept = BMD_pathway,
#       lty = 4,
#       color = "red",
#       linewidth = 0.7
#     ) +
#     theme(panel.grid.major.y = element_line(
#       color = "lightgray",
#       linewidth = 0.5,
#       linetype = 2
#     ),
#     axis.text.y = element_text(size = 6),
#     axis.text.x = element_text(size = 8))
#   
#   return(p)
# }
