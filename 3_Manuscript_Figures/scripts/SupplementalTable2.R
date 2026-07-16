
# supplemental table 2

load(file.path(getwd(), "output", "EUT046", "WrangledInput", "WrangledInputData.RData"))


# get top 15 and make table
top15_LU = LU_norm_HALLMARK_select %>% filter(timepoint == "24h") %>% arrange(medianBMD) %>% head(15) %>% 
  select("pathway_LU" = Pathway.Name, "medianBMD_LU" = medianBMD) %>% mutate(medianBMD_LU = round(medianBMD_LU, 2))
top15_SC = Sciensano_norm_HALLMARK_select %>% filter(timepoint == "24h") %>% arrange(medianBMD) %>% head(15) %>% 
  select("pathway_SC" = Pathway.Name, "medianBMD_SC" = medianBMD) %>% mutate(medianBMD_SC = round(medianBMD_SC, 2))
top15_AU = AU_norm_HALLMARK_select %>% filter(timepoint == "24h") %>% arrange(medianBMD) %>% head(15) %>% 
  select("pathway_AU" = Pathway.Name, "medianBMD_AU" = medianBMD) %>% mutate(medianBMD_AU = round(medianBMD_AU, 2))
top15_BPI = BPI_norm_HALLMARK_select %>% filter(timepoint == "24h") %>% arrange(medianBMD) %>% head(15) %>%
  select("pathway_BPI" = Pathway.Name, "medianBMD_BPI" = medianBMD) %>% mutate(medianBMD_BPI = round(medianBMD_BPI, 2))
top15_GU = GU_norm_HALLMARK_select %>% filter(timepoint == "24h") %>% arrange(medianBMD) %>% head(15) %>% 
  select("pathway_GU" = Pathway.Name, "medianBMD_GU" = medianBMD)%>% mutate(medianBMD_GU = round(medianBMD_GU, 2))


total_suptable2 = bind_cols(
  "top" = seq(1, 15, 1),
  top15_LU,
  top15_SC,
  top15_AU,
  top15_BPI,
  top15_GU
)



library(flextable)
library(officer)

ft = flextable(total_suptable2)
ft = autofit(ft)

doc = read_docx()
doc = body_add_flextable(doc, ft)

print(doc, target = file.path(getwd(), "tables", "supplemental_table2R.docx"))

