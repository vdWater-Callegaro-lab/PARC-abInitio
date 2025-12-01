

tpods = data.table::fread(file.path(getwd(), "output", "EUT046", "tpod_calculations_alltimepoints.txt"))

tpod_bootstrapping = data.table::fread(file.path(getwd(), "output", "EUT046", "tpod_bootstrapping_ci_median.txt"))


# generate table 3

library(tidyverse)
library(stringr)

# join original tpod and bootstrapped results
tpod_joined <- tpods %>% 
  select(partner, timepoint, method, tpod_orig) %>%
  left_join(
    tpod_bootstrapping %>%
      select(partner, timepoint, method, tpod_lower, tpod_upper),
    by = c("partner", "timepoint", "method")
  )

# format tpods: original (lower - upper)
tpod_formatted <- tpod_joined %>%
  mutate(
    tpod_ci = str_glue(
      "{round(tpod_orig, 2)} ({round(tpod_lower, 2)}–{round(tpod_upper, 2)})"
    )
  )

## order methods
method_order <- c("5th percentile",
                  "25th ranked gene",
                  "LCRD",
                  "First mode")


# order timepoints
timepoint_order <- c("4h", "8h", "16h", "24h", "48h", "72h") 

# order partners
partner_order <- c("AristotleU", "BPI", "GhentU", "LeidenU", "Sciensano")

tpod_ordered <- tpod_formatted %>%
  mutate(
    method    = factor(method, levels = method_order),
    timepoint = factor(timepoint, levels = timepoint_order),
    partner = factor(partner, levels = partner_order)
  ) %>%
  arrange(timepoint, method)

# generate final table
tpod_table_final <- tpod_ordered %>%
  select(timepoint, method, partner, tpod_ci) %>%
  pivot_wider(
    names_from  = partner,
    values_from = tpod_ci
  ) %>%
  arrange(timepoint, method)

tpod_table_final


library(flextable)
library(officer)

ft = flextable(tpod_table_final)
ft = autofit(ft)

doc = read_docx()
doc = body_add_flextable(doc, ft)

print(doc, target = file.path(getwd(), "tables", "table3R.docx"))


