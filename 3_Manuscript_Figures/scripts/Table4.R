

tpods = data.table::fread(file.path(getwd(), "output", "EUT046", "tpod_calculations_alltimepoints.txt"))

tpod_bootstrapping = data.table::fread(file.path(getwd(), "output", "EUT046", "tpod_bootstrapping_ci_median.txt"))


# generate table 3

library(tidyverse)
library(stringr)

# join original tpod and bootstrapped results
tpod_joined <- tpods %>% 
  select(analysis_summary, timepoint, method, tpod_orig) %>%
  left_join(
    tpod_bootstrapping %>%
      select(analysis_summary, timepoint, method, tpod_lower, tpod_upper),
    by = c("analysis_summary", "timepoint", "method")
  )

# format tpods: original (lower - upper)
tpod_formatted <- tpod_joined %>%
  mutate(
    tpod_ci = str_glue(
      "{round(tpod_orig, 1)} ({round(tpod_lower, 1)}–{round(tpod_upper, 1)})"
    )
  )

## order methods
method_order <- c("5th percentile",
                  "25th ranked gene",
                  "First mode",
                  "Kneedle",
                  "LCRD")


# order timepoints
timepoint_order <- c("4h", "8h", "16h", "24h", "48h", "72h") 

analysis_summary_order <- c("BMDE-noWTT-CPM-RF-S5", "BMDE-WTT-CPM-RF-S0", "DRO-Quad-UQ-RF-S0", "DRO-Quad-VST-C10-S0", "DRO-Quad-VST-RF-S0")

tpod_ordered <- tpod_formatted %>%
  mutate(
    method    = factor(method, levels = method_order),
    timepoint = factor(timepoint, levels = timepoint_order)
  ) %>%
  arrange(timepoint, method)

tpod_mean <- tpod_joined %>%
  group_by(timepoint, method) %>%
  summarise(mean_tpod = mean(tpod_orig, na.rm = TRUE), .groups = "drop") %>%
  mutate(mean_tpod_fmt = sprintf("%.1f", mean_tpod),
         timepoint = factor(timepoint, levels = timepoint_order))


# generate final table
tpod_table_final <- tpod_ordered %>%
  select(timepoint, method, analysis_summary, tpod_ci) %>%
  pivot_wider(
    names_from  = analysis_summary,
    values_from = tpod_ci
  ) %>%
  left_join(
    tpod_mean %>% select(timepoint, method, mean_tpod_fmt),
    by = c("timepoint", "method")
  ) %>%
  arrange(timepoint, method) %>%
  select(timepoint, method, all_of(analysis_summary_order), mean_tPOD = mean_tpod_fmt)

tpod_table_final


library(flextable)
library(officer)

ft = flextable(tpod_table_final)
ft = autofit(ft)

doc = read_docx()
doc = body_add_flextable(doc, ft)

print(doc, target = file.path(getwd(), "tables", "table4R.docx"))


