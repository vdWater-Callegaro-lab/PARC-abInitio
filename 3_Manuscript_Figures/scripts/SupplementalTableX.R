

# supplemental table X

nselected_combined <- fread("output/EUT046/nselected_crg_after_permutations_combined.txt")


# set timepoint levels
timepoint_levels <- c("4h", "8h", "16h", "24h", "48h", "72h")

table_summary <- nselected_combined %>%
  mutate(timepoint = factor(timepoint, levels = timepoint_levels)) %>%
  group_by(partner, timepoint) %>%
  summarise(
    n_permutations = n_distinct(permutation),
    n_nonzero = sum(n_selected > 0, na.rm = TRUE),
    
    median_positive = if (
      any(n_selected > 0, na.rm = TRUE)
    ) {
      median(n_selected[n_selected > 0], na.rm = TRUE)
    } else {
      NA_real_
    },
    
    min_positive = if (
      any(n_selected > 0, na.rm = TRUE)
    ) {
      min(n_selected[n_selected > 0], na.rm = TRUE)
    } else {
      NA_real_
    },
    
    max_positive = if (
      any(n_selected > 0, na.rm = TRUE)
    ) {
      max(n_selected[n_selected > 0], na.rm = TRUE)
    } else {
      NA_real_
    },
    
    .groups = "drop"
  ) %>%
  mutate(
    selection_frequency = sprintf(
      "%d/%d",
      n_nonzero,
      n_permutations
    ),
    magnitude = case_when(
      n_nonzero == 0 ~ "–",
      TRUE ~ sprintf(
        "%g (%g–%g)",
        median_positive,
        min_positive,
        max_positive
      )
    )
  ) %>%
  select(
    timepoint,
    partner,
    `Nonzero permutations` = selection_frequency,
    `Median (min–max)` = magnitude
  ) %>%
  pivot_longer(
    cols = c(`Nonzero permutations`, `Median (min–max)`),
    names_to = "Statistic",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = partner,
    values_from = Value
  ) %>%
  arrange(
    timepoint,
    factor(
      Statistic,
      levels = c("Nonzero permutations", "Median (min–max)")
    )
  )

table_summary



library(flextable)
library(officer)

ft = flextable(table_summary)
ft = autofit(ft)

doc = read_docx()
doc = body_add_flextable(doc, ft)

print(doc, target = file.path(getwd(), "tables", "supplementalTableXR.docx"))