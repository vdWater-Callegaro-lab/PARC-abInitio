
source("R/render_bmd_report.R")

# Params to override: 
#   params:
#     study_id: EUT046
#     input_dir: input/LeidenU
#     output_dir: output/LeidenU
#     run_bmd_normalized_counts: true
#     run_bmd_log2fc: true
#     run_bmd_pathway_scores: true
#     have_median_pathway_scores: true
#     have_trend_test: true
#     files:
#       value:
#         sample_metadata: sample_metadata.txt
#         normalized_counts: normalized_counts.txt
#         log2fc_table: log2fc_table.txt
#         pathway_overview: pathway_overview_hallmarks.txt
#         bmd_gene_normalized_counts: BMD_output_normalized_counts.txt
#         bmd_gene_log2fc: BMD_output_log2fc.txt
#         bmd_pathway_scores: BMD_output_pathway_scores.txt
#         median_bmd_normalized_counts: medianBMD_pathway_normalized_counts_hallmarks.txt
#         median_bmd_log2fc: medianBMD_pathways_log2fc.txt
#         tt_normalized_counts: WTT_normalized_counts.txt
#         tt_log2fc: WTT_log2fc.txt
#         tt_pathway_scores: WTT_pathway_scores.txt



# render report LeidenU
render_bmd_report(
  rmd_path = "Generalized_BMD_Reporting_Template.Rmd",
  output_basename = paste0(format(Sys.Date(), "%Y%m%d"), "_LeidenU_EUT046_Generalized_BMD_Reporting.txt")
)



# render report Sciensano
render_bmd_report(
  rmd_path = "Generalized_BMD_Reporting_Template.Rmd",
  output_basename = paste0(format(Sys.Date(), "%Y%m%d"), "_Sciensano_EUT046_Generalized_BMD_Reporting.txt"),
  reports_dir = "reports",
  params_override = list(
    study_id = "EUT046",
    input_dir = "input/Sciensano",
    output_dir = "output/Sciensano",
    run_bmd_normalized_counts = TRUE,
    run_bmd_log2fc = FALSE,
    run_bmd_pathway_scores = TRUE,
    have_median_pathway_scores = TRUE,
    have_trend_test = TRUE,
    files = list(
      bmd_gene_normalized_counts = "BMD_output_normalized_counts_incltp.txt",
      bmd_pathway_scores = "BMD_output_NES_final_incltp.txt",
      median_bmd_normalized_counts = "medianBMD_pathways_normalized_counts_finalfiltered_incltp.txt",
      pathway_overview = "pathway_overview.txt",
      tt_normalized_counts = "BMD_output_normalized_counts_incltp.txt",
      tt_pathway_scores = "BMD_output_NES_final_incltp.txt"
    )
  )
)



# render report BPI
render_bmd_report(
  rmd_path = "Generalized_BMD_Reporting_Template.Rmd",
  output_basename = paste0(format(Sys.Date(), "%Y%m%d"), "_BPI_EUT046_Generalized_BMD_Reporting.txt"),
  reports_dir = "reports",
  params_override = list(
    study_id = "EUT046",
    input_dir = "input/BPI",
    output_dir = "output/BPI",
    run_bmd_normalized_counts = TRUE,
    run_bmd_log2fc = FALSE,
    run_bmd_pathway_scores = FALSE,
    have_median_pathway_scores = TRUE,
    have_trend_test = TRUE,
    files = list(
      bmd_gene_normalized_counts = "ALL_FilteredDRomicV3.txt",
      median_bmd_normalized_counts = "All_WPOutV3.txt",
      pathway_overview = "pathway_overview_wikipathways.txt",
      tt_normalized_counts = "All_FilteredDRomicV3.txt",
      normalized_counts = "normalized_counts.txt"
    )
  )
)


# AristotleU
render_bmd_report(
  rmd_path = "Generalized_BMD_Reporting_Template.Rmd",
  output_basename = paste0(format(Sys.Date(), "%Y%m%d"), "_AristotleU_EUT046_Generalized_BMD_Reporting.txt"),
  reports_dir = "reports",
  params_override = list(
    study_id = "EUT046",
    input_dir = "input/AristotleU",
    output_dir = "output/AristotleU",
    run_bmd_normalized_counts = TRUE,
    run_bmd_log2fc = FALSE,
    run_bmd_pathway_scores = FALSE,
    have_median_pathway_scores = TRUE,
    have_trend_test = TRUE,
    files = list(
      bmd_gene_normalized_counts = "EUT046_log2cpm_All_BMD_output_incltp.txt",
      median_bmd_normalized_counts = "EUT046_KEGG_Pathway_Medians_BMD_output_incltp.txt",
      pathway_overview = "EUT046_PathwayIdentifiers.txt",
      tt_normalized_counts = "EUT046_log2cpm_All_BMD_output_incltp.txt",
      normalized_counts = "EUT046_Preprocessed_log2cpm_Universe.txt",
      sample_metadata = "EUT046_metadata_final.txt"
    )
  )
)



# GhentU
render_bmd_report(
  rmd_path = "Generalized_BMD_Reporting_Template.Rmd",
  output_basename = paste0(format(Sys.Date(), "%Y%m%d"), "_GhentU_EUT046_Generalized_BMD_Reporting.txt"),
  reports_dir = "reports",
  params_override = list(
    study_id = "EUT046",
    input_dir = "input/GhentU",
    output_dir = "output/GhentU",
    run_bmd_normalized_counts = TRUE,
    run_bmd_log2fc = FALSE,
    run_bmd_pathway_scores = FALSE,
    have_median_pathway_scores = TRUE,
    have_trend_test = TRUE,
    files = list(
      bmd_gene_normalized_counts = "BMD_output_normalized_counts.txt",
      median_bmd_normalized_counts = "medianBMD_pathways_normalized_counts_filteredonn.txt",
      pathway_overview = "pathway_overview.txt",
      tt_normalized_counts = "BMD_output_normalized_counts.txt"
    )
  )
)

