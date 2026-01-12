

# Generate table 5 (pathway-based modeling vs median BMC aggregation)

# Settings

## Libraries
library(tidyverse)
library(data.table)


## Load data
timepoint_levels = c("4h", "8h", "16h", "24h", "48h", "72h")

load(file.path(getwd(), "output", "EUT046", "WrangledInput", "WrangledInputData.RData")) # = filtered data 


## create table
table5 = tibble(timepoint = timepoint_levels) %>% left_join(
  LU_norm_modules_select %>% filter(Pathway.Name == "hRPTECTERT1_35") %>%
    select(timepoint, "Median BMC (LeidenU)" = medianBMD)
) %>%
  left_join(
    LU_BMD_pathway_select %>% filter(pathway == "hRPTECTERT1_35") %>%
      select(timepoint, "Pathway Score (EGs) Based" = finalBMD)
  ) %>%
  left_join(
    Sciensano_norm_HALLMARK_select %>% filter(Pathway.Name == "HALLMARK_P53_PATHWAY") %>%
      select(timepoint, "Median BMC (Sciensano)" = medianBMD)
  ) %>%
  left_join(
    Sciensano_BMD_pathway_select %>% filter(pathway == "HALLMARK_P53_PATHWAY") %>%
      select(timepoint, "Pathway Score (NES) Based" = finalBMD)
  )



library(flextable)
library(officer)

ft = flextable(table5)
ft = autofit(ft)

doc = read_docx()
doc = body_add_flextable(doc, ft)

print(doc, target = file.path(getwd(), "tables", "table5R.docx"))
