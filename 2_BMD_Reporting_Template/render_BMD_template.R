


# run_BMD_report.R
library(rmarkdown)


rmarkdown::render(
  "scripts/Generalized_BMD_Reporting_Template.Rmd",
  params = list(
    study_id = "EUT046",
    input_dir = "input/EUT046/newAnalysisIQR3",
    output_dir = "output/EUT046",
    figures_dir = "figures/EUT046",
    run_bmd_log2fc = FALSE,
    have_williams_trend_test = TRUE,
    run_bmd_network_scores = FALSE,
    have_uploaded_network_scores = FALSE,
    files = list( # <-- this entire list is the 'value'
      sample_metadata = "sample_metadata.txt",
      normalized_counts = "normalized_counts.txt",
      log2fc_table = "log2fc_table.txt",
      pathway_overview = "pathway_overview.txt",
      bmd_gene_log2cpm = "BMD_output_normalized_counts.txt",
      bmd_gene_log2fc  = "BMD_output_log2fc.txt",
      bmd_pathway_scores = "BMD_output_network_scores.txt",
      pathway_median_bmd_from_genes = "medianBMD_pathways_normalized_counts.txt",
      wtt_log2cpm = "WTT_normalized_counts.txt",
      wtt_log2fc  = "WTT_log2fc.txt"
    )
    # columns = list( ... )  # same idea if you want to override the mappings
  )
)
