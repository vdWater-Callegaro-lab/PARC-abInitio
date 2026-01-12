
# Table 4
## Include obtained BMCs for HALLMARK p53 & p53 pathways of choice

### Load packages
library(tidyverse)
library(flextable)
library(officer)

### Load data
load(file.path(getwd(), "output", "EUT046", "WrangledInput", "WrangledInputData.RData")) # = filtered data 




#### HALLMARK
partner_order = c("AU", "BPI", "GU", "LU", "SC")
timepoint_order = c("4h", "8h", "16h", "24h", "48h", "72h")


result_HALLMARK_p53 <- norm_HALLMARK_combined %>%
  filter(Pathway.Name == "HALLMARK_P53_PATHWAY") %>%
  select(Pathway.Name, timepoint, medianBMD, partner) %>%
  # enforce factor levels for ordering
  mutate(
    partner   = factor(partner, levels = partner_order),
    timepoint = factor(timepoint, levels = timepoint_order),
    medianBMD = round(medianBMD, digits = 3)
  ) %>%
  # ensure all partner–timepoint combinations exist
  complete(timepoint, partner) %>%
  # reshape so partners become columns
  pivot_wider(
    names_from  = partner,
    values_from = medianBMD
  ) %>%
  # order rows explicitly
  arrange(timepoint) %>%
  select(Pathway.Name, everything()) %>%
  filter(!is.na(Pathway.Name)) %>%
  dplyr::rename(
    "AristotleU" = AU,
    "GhentU" = GU,
    "LeidenU" = LU,
    "Sciensano" = SC
  )


# save

ft = flextable(result_HALLMARK_p53)
ft = autofit(ft)

doc = read_docx()
doc = body_add_flextable(doc, ft)

print(doc, target = file.path(getwd(), "tables", "table4aR.docx"))



### Pathway of choice
p53_pathways = c("p53 signaling pathway", "p53 transcriptional gene network", "HALLMARK_P53_PATHWAY", "hRPTECTERT1_35", "GENOMARK84", "TGxDDI64")

result_p53_choice = norm_pathwaysChoice_combined %>%
  filter(Pathway.Name %in% p53_pathways) %>%
  select(Pathway.Name, timepoint, medianBMD) %>%
  mutate(
    timepoint = factor(timepoint, levels = timepoint_order),
    Pathway.Name = factor(Pathway.Name, levels = p53_pathways),
    medianBMD = round(medianBMD, digits = 3)
  ) %>%
  complete(timepoint, Pathway.Name) %>%
  pivot_wider(
    names_from = Pathway.Name, 
    values_from = medianBMD
  ) %>%
  arrange(timepoint)





ft = flextable(result_p53_choice)
ft = autofit(ft)

doc = read_docx()
doc = body_add_flextable(doc, ft)

print(doc, target = file.path(getwd(), "tables", "table4bR.docx"))



