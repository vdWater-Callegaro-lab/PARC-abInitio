
# Generate table with Dose Responsive Genes + Retained genes after Post Model Filters

# Settings

## Libraries
library(tidyverse)
library(data.table)


## Directories
studyNR = "EUT046"
inputDir = file.path(here(), "input", studyNR)
outputDir = file.path(here(), "output", studyNR)


## Load data

timepoint_levels = c("4h", "8h", "16h", "24h", "48h", "72h")

### Prefilter output
AU_pref = fread(file.path(inputDir, "Aristotle", paste0(studyNR, "_log2cpm_All_BMD_output.txt")), check.names = TRUE) %>% 
  mutate(timepoint = factor(
  str_extract(Analysis, "\\d{1,2}h"),
  levels = timepoint_levels)
  ) %>%
  filter(adjpvalue < 0.01)

BPI_pref = fread(file.path(inputDir, "BPI", "All_UnfilteredDRomicV3.txt"), check.names = TRUE) %>%
  filter(adjpvalue < 0.01) %>%
  mutate(timepoint = factor(timepoint, levels = timepoint_levels)) 


folder_path_GU <- "input/EUT046/Ghent/Ab_initio_UGent/Output_046/Bootstrap_unfiltered"
GU_pref <- list.files(folder_path_GU, pattern = "\\.txt$", full.names = TRUE) %>%
  map_dfr(~ {
    fread(.x, check.names = TRUE) %>%
      mutate(
        timepoint = factor(paste0(
          stringr::str_extract(basename(.x), "(?<=Time_)\\d+"),
          "h"
        ), levels = timepoint_levels
      )
      )
  }) %>%
  filter(adjpvalue < 0.01)


LU_pref = fread(file.path(inputDir, "Leiden", paste0(studyNR, "_log2cpm_WTT_output.txt")), check.names = TRUE) %>% 
  mutate(timepoint = factor(
  str_extract(Analysis, "\\d{1,2}h"), 
  levels = timepoint_levels)
  ) %>% filter(Adjusted.P.Value < 0.05) %>%
  filter(Probe.ID != "NA_NA")


SC_pref = fread(file.path(inputDir, "Sciensano", "BMD_output_normalized_counts.txt"), check.names = TRUE) %>% 
  mutate(timepoint = factor(
  str_extract(Analysis, "\\d{1,2}h"), 
  levels = timepoint_levels)
)



### After post model filters
load(file.path(getwd(), "output", "EUT046", "WrangledInput", "WrangledInputData.RData")) # = filtered data 




# dose responsive genes
drg_afterprefilter = tibble(
  "DRomics_CPM_QT" = AU_pref %>% pull(timepoint) %>% table() %>% as.numeric(),
  "DRomics_log2Internal_QT" = BPI_pref %>% pull(timepoint) %>% table() %>% as.numeric(),
  "DRomics_VST_QT" = GU_pref %>% pull(timepoint) %>% table() %>% as.numeric(),
  "BMDExpress_log2CPM_noWTT" = LU_pref %>% pull(timepoint) %>% table() %>% as.numeric(),
  "BMDExpress_log2CPM_WTT" = SC_pref %>% pull(timepoint) %>% table() %>% as.numeric()
) %>%
  rowwise() %>%
  mutate(
    mean_cv = sprintf(
      "%.2f (%.1f%%)",
      mean(c_across(1:5), na.rm = TRUE),
      sd(c_across(1:5), na.rm = TRUE) /
        mean(c_across(1:5), na.rm = TRUE) * 100
    )
  ) %>%
  ungroup() %>%
  mutate(timepoint = factor(timepoint_levels))

drg_afterprefilter

data.table::fwrite(drg_afterprefilter, file.path(getwd(), "output", "EUT046", "drg_afterprefilter.txt"))




# retained after post model filters
retained_pmf = tibble(
  "DRomics_CPM_QT" = AU_norm_BMD_select %>% pull(timepoint) %>% table() %>% as.numeric(),
  "DRomics_log2Internal_QT" = BPI_norm_BMD_select %>% pull(timepoint) %>% table() %>% as.numeric(),
  "DRomics_VST_QT" = GU_norm_BMD_select %>% pull(timepoint) %>% table() %>% as.numeric(),
  "BMDExpress_log2CPM_noWTT" = LU_norm_BMD_select %>% pull(timepoint) %>% table() %>% as.numeric(),
  "BMDExpress_log2CPM_WTT" = Sciensano_norm_BMD_select %>% pull(timepoint) %>% table() %>% as.numeric()
) %>%
  rowwise() %>%
  mutate(
    mean_cv = sprintf(
      "%.2f (%.1f%%)",
      mean(c_across(1:5), na.rm = TRUE),
      sd(c_across(1:5), na.rm = TRUE) /
        mean(c_across(1:5), na.rm = TRUE) * 100
    )
  ) %>%
  ungroup() %>%
  mutate(timepoint = factor(timepoint_levels))



retained_pmf

data.table::fwrite(retained_pmf, file.path(getwd(), "output", "EUT046", "retained_after_postmodelfilters.txt"))




# combine
DRG_PMF_table2 = bind_rows(
  drg_afterprefilter,
  retained_pmf
) %>%
  select(timepoint, everything())




library(flextable)
library(officer)

ft = flextable(DRG_PMF_table2)
ft = autofit(ft)

doc = read_docx()
doc = body_add_flextable(doc, ft)

print(doc, target = file.path(getwd(), "tables", "table2R.docx"))



#  create long table
drg_afterprefilter_long <- drg_afterprefilter %>%
  select(-mean_cv) %>%
  pivot_longer(cols = -timepoint, names_to = "analysis_summary", values_to = "value") %>%
  mutate(step = "DRG")


retained_pmf_long <- retained_pmf %>%
  select(-mean_cv) %>%
  pivot_longer(cols = -timepoint, names_to = "analysis_summary", values_to = "value") %>%
  mutate(step = "PMFG")


drg_pmf_combined_long = bind_rows(drg_afterprefilter_long, retained_pmf_long)
data.table::fwrite(drg_pmf_combined_long, file.path(outputDir, "DRG_PMFG_pertimepoint.txt"))
