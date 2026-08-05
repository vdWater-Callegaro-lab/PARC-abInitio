

# supplemental table X 

nselected_combined <- fread("output/EUT046/nselected_crg_after_permutations_combined.txt")
drgs <- fread("output/EUT046/drg_afterprefilter.txt")


# set timepoint levels
timepoint_levels <- c("4h", "8h", "16h", "24h", "48h", "72h")



# perm_results should contain:
# timepoint, n_selected, permutation, analysis
#
# drgs should contain:
# timepoint + one column per analysis, with values = observed CRG counts

# 1. Convert observed CRG counts from wide to long format
observed_results <- drgs %>%
  select(-mean_cv) %>%
  pivot_longer(
    cols = -timepoint,
    names_to = "analysis",
    values_to = "observed_n_selected"
  ) %>%
  mutate(
    timepoint = as.character(timepoint),
    analysis = as.character(analysis),
    observed_n_selected = as.numeric(observed_n_selected)
  )

# 2. Standardize permutation results
perm_results_clean <- nselected_combined %>%
  mutate(
    timepoint = as.character(timepoint),
    analysis = as.character(analysis),
    n_selected = as.numeric(n_selected),
    nonzero = n_selected > 0
  )

# Optional but useful: check whether all analysis/time point combinations match
missing_observed <- perm_results_clean %>%
  distinct(analysis, timepoint) %>%
  anti_join(
    observed_results %>% distinct(analysis, timepoint),
    by = c("analysis", "timepoint")
  )

if (nrow(missing_observed) > 0) {
  print(missing_observed)
  stop("Some analysis/timepoint combinations in perm_results are missing from drgs.")
}

# 3. Generate supplemental table with empirical p-value
supp_table_empirical_null <- perm_results_clean %>%
  left_join(
    observed_results,
    by = c("analysis", "timepoint")
  ) %>%
  group_by(analysis, timepoint, observed_n_selected) %>%
  summarise(
    n_permutations = sum(!is.na(n_selected)),
    n_nonzero_permutations = sum(nonzero, na.rm = TRUE),
    percent_nonzero_permutations = 100 * mean(nonzero, na.rm = TRUE),
    
    mean_null_crgs = mean(n_selected, na.rm = TRUE),
    median_null_crgs = median(n_selected, na.rm = TRUE),
    q25_null_crgs = quantile(n_selected, 0.25, na.rm = TRUE),
    q75_null_crgs = quantile(n_selected, 0.75, na.rm = TRUE),
    min_null_crgs = min(n_selected, na.rm = TRUE),
    max_null_crgs = max(n_selected, na.rm = TRUE),
    p95_null_crgs = quantile(n_selected, 0.95, na.rm = TRUE),
    
    n_permutations_ge_observed = sum(n_selected >= observed_n_selected, na.rm = TRUE),
    empirical_p_value = (1 + n_permutations_ge_observed) / (1 + n_permutations),
    
    .groups = "drop"
  ) %>%
  mutate(
    percent_nonzero_permutations = round(percent_nonzero_permutations, 1),
    mean_null_crgs = round(mean_null_crgs, 1),
    median_null_crgs = round(median_null_crgs, 1),
    q25_null_crgs = round(q25_null_crgs, 1),
    q75_null_crgs = round(q75_null_crgs, 1),
    p95_null_crgs = round(p95_null_crgs, 1),
    null_iqr_crgs = paste0(q25_null_crgs, "–", q75_null_crgs),
    null_range_crgs = paste0(min_null_crgs, "–", max_null_crgs),
    empirical_p_value = signif(empirical_p_value, 3)
  ) %>%
  select(
    analysis,
    timepoint,
    observed_n_selected,
    percent_nonzero_permutations,
    median_null_crgs,
    null_range_crgs,
    empirical_p_value
  ) %>%
  arrange(
    analysis,
    as.numeric(str_remove_all(timepoint, "[^0-9.]"))
  ) %>%
  dplyr::rename(
    `Workflow` = analysis,
    `Time point` = timepoint,
    `Observed CRGs` = observed_n_selected,
    `Permutations with ≥1 CRG (%)` = percent_nonzero_permutations,
    `Median null CRGs` = median_null_crgs,
    `Null CRG range` = null_range_crgs,
    `Empirical p-value` = empirical_p_value
  )


supp_table_empirical_null



library(flextable)
library(officer)

ft = flextable(supp_table_empirical_null)
ft = autofit(ft)

doc = read_docx()
doc = body_add_flextable(doc, ft)

print(doc, target = file.path(getwd(), "tables", "supplementalTableXR.docx"))

