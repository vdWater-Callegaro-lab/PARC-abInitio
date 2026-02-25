

# BOOTSTRAPPING TO IDENTIFY CONFIDENCE INTERVAL IN TPOD CALCULATIONS

## LOAD PACKAGES
library(tidyverse)
library(purrr)

## LOAD FUNCTION TO CALUCATE TPODS
source(file.path(getwd(), "functions.R"))


## LOAD DATA
load(file.path(getwd(), "output", "EUT046", "WrangledInput", "WrangledInputData.RData"))



## ORIGINAL TPODS
tpod_orig_all = data.table::fread(file.path(getwd(), "output", "EUT046", "tpod_calculations_alltimepoints.txt"))




## BOOTSTRAPPING

### GENERATE BOOTSTRAP FUNCTION
bootstrap_tpods <- function(df_bmd,
                            n_boot = 1000,
                            seed = 123) {
  set.seed(seed)
  
  map_dfr(seq_len(n_boot), function(b) {
    
    # 1) Resample genes WITHIN each timepoint
    df_boot <- df_bmd %>%
      group_by(timepoint) %>%
      group_modify(~ {
        if (nrow(.x) < 5) {
          # if too few genes, return empty tibble for this timepoint
          .x[0, ]
        } else {
          dplyr::slice_sample(.x, n = nrow(.x), replace = TRUE)
          # for older dplyr: dplyr::sample_n(.x, size = nrow(.x), replace = TRUE)
        }
      }) %>%
      ungroup()
    
    # 2) Run your existing tPOD functions on the bootstrapped df
    res_5   <- get_percentile_bmd(df_boot, prob = 0.05)        %>% mutate(method = "p5")
    res_25  <- get_nth_ranked_bmd(df_boot, n = 25)       %>% mutate(method = "rank25")
    res_lcrd <- get_LCRD_bmd(df_boot)      %>% mutate(method = "LCRD")
    res_mode <- get_first_mode_bmd(df_boot)%>% mutate(method = "first_mode")
    res_kneedle <- get_kneedle_bmd(df_boot) %>% mutate(method = "kneedle")
    
    # 3) Stack all methods and tag with bootstrap iteration
    bind_rows(res_5, res_25, res_lcrd, res_mode, res_kneedle) %>%
      mutate(boot = b, .before = 1)
  })
}



### RUN BOOTSTRAP FUNCTION
bootstrapped_tpods_LU = bootstrap_tpods(df_bmd = LU_norm_BMD_select %>% 
                                          dplyr::rename("bmd" = finalBMD), n_boot = 1000, seed = 123)

bootstrapped_tpods_SC = bootstrap_tpods(df_bmd = Sciensano_norm_BMD_select %>% 
                                          dplyr::rename("bmd" = finalBMD), n_boot = 1000, seed = 123)

bootstrapped_tpods_BPI = bootstrap_tpods(df_bmd = BPI_norm_BMD_select %>% 
                                          dplyr::rename("bmd" = finalBMD), n_boot = 1000, seed = 123)

bootstrapped_tpods_GU = bootstrap_tpods(df_bmd = GU_norm_BMD_select %>% 
                                          dplyr::rename("bmd" = finalBMD), n_boot = 1000, seed = 123)

bootstrapped_tpods_AU = bootstrap_tpods(df_bmd = AU_norm_BMD_select %>% 
                                          dplyr::rename("bmd" = finalBMD), n_boot = 1000, seed = 123)



### COMBINE RESULTS
boot_all <- list(
  `BMDE-noWTT-CPM-RF-S5`    = bootstrapped_tpods_LU,
  `BMDE-WTT-CPM-RF-S0`  = bootstrapped_tpods_SC,
  `DRO-Quad-VST-RF-S0`     = bootstrapped_tpods_GU,
  `DRO-Quad-VST-C10-S0`        = bootstrapped_tpods_BPI,
  `DRO-Quad-UQ-RF-S0` = bootstrapped_tpods_AU
) %>%
  imap_dfr(~ mutate(.x, analysis_summary = .y))   # add analysis_summary column from list names

data.table::fwrite(boot_all, file.path(getwd(), "output", "EUT046", "tpod_bootstrapping_results.txt"), sep = "\t")





### GET CONFIDENCE INTERVALS PER METHOD X PARTNER X TIMEPOINT
tpod_ci_all <- boot_all %>%
  group_by(analysis_summary, timepoint, method) %>%
  summarise(
    tpod_median = median(tpod, na.rm = TRUE),
    tpod_lower  = quantile(tpod, 0.025, na.rm = TRUE),
    tpod_upper  = quantile(tpod, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    method = factor(method,
                    levels = c("p5","rank25","LCRD","first_mode", "kneedle"),
                    labels = c("5th percentile","25th ranked gene","LCRD","First mode", "Kneedle"))
  )


data.table::fwrite(tpod_ci_all, file.path(getwd(), "output", "EUT046", "tpod_bootstrapping_ci_median.txt"), sep = "\t")


