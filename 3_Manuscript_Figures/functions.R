
# functions


# BMD accumulation plot
BMDaccumulationPlot <- function(BMDoutput, partner) {
  BMDoutput %>%
    group_by(timepoint) %>%
    arrange(finalBMD) %>%
    mutate(ranked_gene = row_number()) %>%
    ggplot(aes(finalBMD, y = ranked_gene)) +
    geom_point(stat = "identity", aes(color = timepoint, shape = timepoint), size = 0.8) +
    labs(title = partner,
         y = "Accumulation",
         x = "BMC") +
    scale_color_manual(values = c("#E6E6E6", "#D1D1D1","#BBBBBB","#A0A0A0","#7F7F7F","#4D4D4D"))
}




# plot BMD accumulation genes zoomed in
BMDaccumulationPlot_zoom <- function(BMDoutput, partner) {
  BMDoutput %>%
    group_by(timepoint) %>%
    arrange(finalBMD) %>%
    mutate(ranked_gene = row_number()) %>%
    ggplot(aes(finalBMD, y = ranked_gene)) +
    geom_point(stat = "identity", aes(color = timepoint, shape = timepoint), size = 1) +
    labs(y = "Accumulation", 
         x = "BMC") +
    ylim(c(0, 250)) + 
    xlim(c(0, 2)) +
    scale_color_manual(values = c("#E6E6E6", "#D1D1D1","#BBBBBB","#A0A0A0","#7F7F7F","#4D4D4D"))
}




BMDaccumulationPlot_zoom_HL5th <- function(BMDoutput, HLgene4hr, HLgene8hr, HLgene16hr, HLgene24hr, HLgene48hr, HLgene72hr){
    plotdata = BMDoutput %>% 
      separate(gene_symbol_entrez_id, into = c("gene_symbol", "entrez_id"), sep = "_") %>%
      group_by(timepoint) %>%
      arrange(finalBMD) %>%
      mutate(ranked_gene = row_number())
    
    highlight_data = plotdata %>% filter(gene_symbol == HLgene4hr & timepoint == "4hr" |
                                           gene_symbol == HLgene8hr & timepoint == "8hr" |
                                           gene_symbol == HLgene16hr & timepoint == "16hr" |
                                           gene_symbol == HLgene24hr & timepoint == "24hr" |
                                           gene_symbol == HLgene48hr & timepoint == "48hr" |
                                           gene_symbol == HLgene72hr & timepoint == "72hr")
    
    ggplot(plotdata, aes(finalBMD, y = ranked_gene)) +
      geom_point(stat = "identity", aes(color = timepoint, shape = timepoint), size = 1) +
      geom_point(data = highlight_data, aes(x = finalBMD, y=ranked_gene, shape = timepoint), color = "red", size = 2) +
      labs(y = "Accumulation", 
           x = "BMD") +
      ylim(c(0, 250)) + 
      xlim(c(0, 2)) +
      theme(legend.position = "none") +
      scale_color_manual(values = c("#E6E6E6", "#D1D1D1","#BBBBBB","#A0A0A0","#7F7F7F","#4D4D4D"))
}




# plot median BMD accumulation of pathways
BMDaccumulationPlot_pathway <- function(BMDoutput, partner) {
  BMDoutput %>%
    group_by(timepoint) %>%
    arrange(medianBMD) %>%
    mutate(ranked_pathway = row_number()) %>%
    ggplot(aes(medianBMD, y = ranked_pathway)) +
    geom_point(stat = "identity", aes(color = timepoint, shape = timepoint), size = 1.5) +
    labs(title = partner,
         y = "Accumulation", 
         x = "BMC") +
    scale_color_manual(values = c("#E6E6E6", "#D1D1D1","#BBBBBB","#A0A0A0","#7F7F7F","#4D4D4D"),
                       breaks = c("4hr", "8hr", "16hr", "24hr", "48hr", "72hr")) +
    scale_shape_manual(values = c(16, 17, 15, 3, 7, 8),
                       breaks = c("4hr", "8hr", "16hr", "24hr", "48hr", "72hr"))
}




# plot median BMD accumulation of pathways + highlight pathway_oi
BMDaccumulationPlot_pathway_highlight <- function(BMDoutput, partner, pathway_oi) {
  plot_data =  BMDoutput %>%
    group_by(timepoint) %>%
    arrange(medianBMD) %>%
    mutate(ranked_pathway = row_number()) 
  
  highlight_data = plot_data %>% filter(Pathway.Name %in% pathway_oi)

    ggplot(plot_data, aes(x= medianBMD, y = ranked_pathway)) +
    geom_point(stat = "identity", aes(color = timepoint, shape = timepoint), size = 1.5) +
    geom_point(data = highlight_data, aes(x = medianBMD, y = ranked_pathway, shape = timepoint), color = "red", size = 3) +
    # geom_text(data = highlight_data, aes(x = medianBMD, y = ranked_pathway, label = Pathway.Name), vjust = -1, hjust = 1, color = "red") +
    labs(title = partner,
         y = "Accumulation",
         x = "BMC") +
    scale_color_manual(values = c("#E6E6E6", "#D1D1D1","#BBBBBB","#A0A0A0","#7F7F7F","#4D4D4D"))
    
}




# plot median BMD accumulation of pathways for all partners combined in one plot + highlight pathway_oi
BMDaccumulationPlot_pathway_highlight_combined <- function(BMDoutput, timepoint_oi, pathway_oi) {
  plot_data = BMDoutput %>%
    filter(timepoint == timepoint_oi) %>%
    group_by(partner) %>%
    arrange(medianBMD) %>%
    mutate(ranked_pathway = row_number()) 
  
  highlight_data = plot_data %>% filter(Pathway.Name %in% pathway_oi)
  
  ggplot(plot_data, aes(x= medianBMD, y = ranked_pathway)) +
    geom_point(stat = "identity", aes(color = partner), size = 1.5) +
    geom_point(data = highlight_data, aes(x = medianBMD, y = ranked_pathway, shape = partner), color = "black", size = 3) +
    labs(y = "Accumulation",
         x = "BMC",
         subtitle = timepoint_oi) +
    scale_color_manual(values = c("firebrick", "gray30", "gold1", "blue2", "forestgreen"),
                       breaks = c("AU", "BPI", "GU", "LU", "SC"),
                       labels = c("Aristotle University (AU)", "BPI", "Ghent University (GU)", "Leiden University (LU)", "Sciensano (SC)")) +
    scale_shape_manual(values = c(0, 8, 1, 2, 6)) +
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
BMDaccumulationPlot_pwscores_highlight <- function(BMDoutput, HLpw, partner){
  plotdata = BMDoutput %>%
    group_by(timepoint) %>%
    arrange(finalBMD) %>%
    mutate(ranked_pathway = row_number())
  
  highlight_data = plotdata %>% filter(pathway == HLpw)
  
  ggplot(plotdata, aes(finalBMD, y = ranked_pathway)) + 
    geom_point(stat = "identity", aes(color = timepoint, shape = timepoint), size = 1) +
    geom_point(data = highlight_data, aes(x = finalBMD, y=ranked_pathway, shape = timepoint), color = "red", size = 3) +
    labs(title = partner,
         y = "Accumulation", 
         x = "BMC") +
    scale_color_manual(values = c("#E6E6E6", "#D1D1D1","#BBBBBB","#A0A0A0","#7F7F7F","#4D4D4D"))

}





# plot dose response curves modules > now does not work anymore; can probably be removed as we do not have this in our figure anymore
BMD_output_plots_pathways = function(BMDoutput, time_oi, BMD_pathways_oi, txgmapr_egs) {
  
  if((BMDoutput %>% filter(timepoint == time_oi) %>% nrow()) == 0){
    print("There are no pathways showing a dose-response")
  }else if((BMDoutput %>% filter(timepoint == time_oi) %>% nrow()) < 5){
    print("NOTE: There are less than 5 pathways showing a dose-response")
  } 
  
  for (pathway in BMD_pathways_oi) {
    BMD_pathway = pathway
    
    BMD_data <-
      BMDoutput %>% filter(patwhay == BMD_pathway &
                                     timepoint == time_oi)
    
    txg_egs_meta_data = txgmapr_egs %>% as.data.frame() %>% # select module of interest from txg data
      filter(module == BMD_module) %>%
      separate(
        sample_id,
        into = c(
          "upload",
          "exposure_type",
          "replicate", 
          "timepoint", 
          "dose"
        ),
        sep="_",
        remove=FALSE
      ) %>%
      mutate(
      #   dose = case_when(
      #   dose_level == "1" ~ 0.1,
      #   dose_level == "2" ~ 0.5,
      #   dose_level == "3" ~ 1,
      #   dose_level == "4" ~ 2.5,
      #   dose_level == "5" ~ 5,
      #   dose_level == "6" ~ 10,
      #   dose_level == "7" ~ 20,
      #   dose_level == "8" ~ 30,
      #   dose_level == "9" ~ 50
      # ),
      timepoint = paste0(timepoint, "hr"),
      dose = as.numeric(dose),
      replicate = gsub("CISPLATIN", "", replicate),
      replicate = paste0("R", replicate)) %>%
      filter(timepoint == time_oi) 
    
    
    BMD_padj = BMD_data$Prefilter.Adjusted.P.Value
    BMD_padj_gene = BMD_data$Min.p.adjust
    # BMD_corEG = BMD_data$corEG
    BMDL = BMD_data$finalBMDL
    BMDU = BMD_data$finalBMDU
    # Cmax_low = Cmax_data$Cmax_low[Cmax_data$COMPOUND == BMD_compound]
    # Cmax_high = Cmax_data$Cmax_high[Cmax_data$COMPOUND == BMD_compound]
    BMD_y = BMD_data$finalBMD
    BMD_model = BMD_data$Best.Model
    BMD_module = BMD_data$module
    
    # if (BMD_data$Max.Fold.Change > 0) {
    #   BMD_logFC = log2(BMD_data$Max.Fold.Change.Absolute.Value)
    # } else {
    #   BMD_logFC = log2(1 / BMD_data$Max.Fold.Change.Absolute.Value)
    # }
    
    HILL_INTERCEPT = BMD_data$Hill.Parameter.Intercept
    HILL_V = BMD_data$Hill.Parameter.v
    HILL_N = BMD_data$Hill.Parameter.n
    HILL_K = BMD_data$Hill.Parameter.k
    
    POWER_CONTROL = BMD_data$Power.Parameter.control
    POWER_SLOPE = BMD_data$Power.Parameter.slope
    POWER_POWER = BMD_data$Power.Parameter.power
    
    LINEAR_0 = BMD_data$Linear.Parameter.beta_0
    LINEAR_1 = BMD_data$Linear.Parameter.beta_1
    
    POLY2_0 = BMD_data$Poly.2.Parameter.beta_0
    POLY2_1 = BMD_data$Poly.2.Parameter.beta_1
    POLY2_2 = BMD_data$Poly.2.Parameter.beta_2
    
    EXP2_A = BMD_data$Exp.2.Parameter.a
    EXP2_B = BMD_data$Exp.2.Parameter.b
    EXP2_SIGN = BMD_data$Exp.2.Parameter.sign
    
    EXP3_A = BMD_data$Exp.3.Parameter.a
    EXP3_B = BMD_data$Exp.3.Parameter.b
    EXP3_D = BMD_data$Exp.3.Parameter.d
    EXP3_SIGN = BMD_data$Exp.3.Parameter.sign
    
    EXP4_A = BMD_data$Exp.4.Parameter.a
    EXP4_B = BMD_data$Exp.4.Parameter.b
    EXP4_C = BMD_data$Exp.4.Parameter.c
    
    EXP5_A = BMD_data$Exp.5.Parameter.a
    EXP5_B = BMD_data$Exp.5.Parameter.b
    EXP5_C = BMD_data$Exp.5.Parameter.c
    EXP5_D = BMD_data$Exp.5.Parameter.d
    
    
    
    HILL   = function(DOSE)
      HILL_INTERCEPT + HILL_V * DOSE ^ HILL_N / (HILL_K ^ HILL_N + DOSE ^ HILL_N)
    POWER  = function(DOSE)
      POWER_CONTROL + POWER_SLOPE * DOSE ^ POWER_POWER
    LINEAR = function(DOSE)
      LINEAR_0 + LINEAR_1 * DOSE
    POLY2  = function(DOSE)
      POLY2_0 + POLY2_1 * DOSE + POLY2_2 * DOSE ^ 2
    EXP2   = function(DOSE)
      EXP2_A * exp(EXP2_SIGN * EXP2_B * DOSE)
    EXP3   = function(DOSE)
      EXP3_A * exp(EXP3_SIGN * (EXP3_B * DOSE) ^ EXP3_D)
    EXP4   = function(DOSE)
      EXP4_A * (EXP4_C - (EXP4_C - 1) * exp(-EXP4_B * DOSE))
    EXP5   = function(DOSE)
      EXP5_A * (EXP5_C - (EXP5_C - 1) * exp(-(EXP5_B * DOSE) ^ EXP5_D))
    
    if (BMD_model == "Hill") {
      BMD_function = HILL
    } else if (BMD_model == "Power") {
      BMD_function = POWER
    } else if (BMD_model == "Linear") {
      BMD_function = LINEAR
    } else if (BMD_model == "Poly 2") {
      BMD_function = POLY2
    } else if (BMD_model == "Exp 2") {
      BMD_function = EXP2
    } else if (BMD_model == "Exp 3") {
      BMD_function = EXP3
    } else if (BMD_model == "Exp 4") {
      BMD_function = EXP4
    } else if (BMD_model == "Exp 5") {
      BMD_function = EXP5
    }
    
    yaxis_lim = c(min(txg_egs_meta_data$eg_score),
                  max(txg_egs_meta_data$eg_score))
    
    xaxis_lim = c(min(txg_egs_meta_data$dose),
                  max(txg_egs_meta_data$dose))
    
    
    p = ggplot(data = txg_egs_meta_data, mapping = aes(x = dose, y = eg_score)) +
      # theme_standard() +
      ggtitle(
        paste0(time_oi, " - ", BMD_module, " - ", BMD_model),
        subtitle = paste0("BMD: ", round(BMD_y, 2), " µM, WTT p-adj: ", round(BMD_padj, 5))
      ) +
      geom_point(aes(group = replicate, shape = replicate), size = 3) +
      stat_function(fun = BMD_function,
                    colour = "blue",
                    size = 1) +
      geom_rect(
        aes(
          xmin = BMDL,
          xmax = BMDU,
          ymin = -Inf,
          ymax = Inf
        ),
        fill = "pink",
        alpha = 0.01
      ) +
      ylim(yaxis_lim) +
      xlab("Concentration (µM)") + ylab("Module eigengene score") +
      scale_shape_manual(name = "Replicate",
                         labels = c("R1", "R2", "R3"),
                         values = rep(seq(21, 25), 5)) +
      geom_vline(xintercept = BMD_y, lty=4, color="blue", linewidth=1) +
      scale_x_log10() 
    
    return(p)
    
  }

}



plot_BMD_genes_moduleBMD = function(BMDmodules, time_oi, module_oi, BMDoutput_genes, pathway_overview) {
  
  BMDL_pathway = BMDmodules %>% filter(timepoint == time_oi & module == module_oi) %>% pull(finalBMDL) %>% unique()
  BMD_pathway = BMDmodules %>% filter(timepoint == time_oi & module == module_oi) %>% pull(finalBMD) %>% unique()
  BMDU_pathway = BMDmodules %>% filter(timepoint == time_oi & module == module_oi) %>% pull(finalBMDU) %>% unique()
  
  genes_in_pathway = pathway_overview %>% filter(module == module_oi) %>% 
    separate(gene_symbol_entrez_id, into = c("gene_symbol", "entrez_id"), sep = "_") %>% 
    pull(gene_symbol) %>% unique()
  
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
