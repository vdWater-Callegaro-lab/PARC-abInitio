

# supplemental table X

nselected_combined <- fread("output/EUT046/nselected_crg_after_permutations_combined.txt")



table_frequency <- nselected_combined %>%
  group_by(partner, timepoint) %>%
  summarise(
    n_permutations = n_distinct(permutation),
    n_nonzero = sum(n_selected > 0, na.rm = TRUE),
    pct_nonzero = 100 * mean(n_selected > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    selection_frequency = sprintf(
      "%d/%d",
      n_nonzero,
      n_permutations
    )
  ) %>%
  select(timepoint, partner, selection_frequency) %>%
  pivot_wider(
    names_from = partner,
    values_from = selection_frequency
  ) %>%
  arrange(timepoint)

table_frequency