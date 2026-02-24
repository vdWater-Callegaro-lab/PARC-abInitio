

library(tidyverse)

source("functions.R")


load(file.path(getwd(), "output", "EUT046", "WrangledInput", "WrangledInputData.RData"))


## GENERATE TABLE 3

# calcualte original tPODs
tpod_orig_LU <- bind_rows(
  get_percentile_bmd(LU_norm_BMD_select %>%
                       dplyr::rename("bmd" = finalBMD), prob = 0.05)        %>% mutate(method = "p5"),
  get_nth_ranked_bmd(LU_norm_BMD_select %>%
                       dplyr::rename("bmd" = finalBMD), n = 25)       %>% mutate(method = "rank25"),
  get_LCRD_bmd(LU_norm_BMD_select %>%
                 dplyr::rename("bmd" = finalBMD))      %>% mutate(method = "LCRD"),
  get_first_mode_bmd(LU_norm_BMD_select %>%
                       dplyr::rename("bmd" = finalBMD)) %>% mutate(method = "first_mode"),
  get_kneedle_bmd(LU_norm_BMD_select %>%
                    dplyr::rename("bmd" = finalBMD)) %>% mutate(method = "kneedle")
)


tpod_orig_SC <- bind_rows(
  get_percentile_bmd(
    Sciensano_norm_BMD_select %>%
      dplyr::rename("bmd" = finalBMD),
    prob = 0.05
  )        %>% mutate(method = "p5"),
  get_nth_ranked_bmd(
    Sciensano_norm_BMD_select %>%
      dplyr::rename("bmd" = finalBMD),
    n = 25
  )       %>% mutate(method = "rank25"),
  get_LCRD_bmd(Sciensano_norm_BMD_select %>%
                 dplyr::rename("bmd" = finalBMD))      %>% mutate(method = "LCRD"),
  get_first_mode_bmd(Sciensano_norm_BMD_select %>%
                       dplyr::rename("bmd" = finalBMD)) %>% mutate(method = "first_mode"),
  get_kneedle_bmd(Sciensano_norm_BMD_select %>%
                    dplyr::rename("bmd" = finalBMD)) %>% mutate(method = "kneedle")
)



tpod_orig_BPI <- bind_rows(
  get_percentile_bmd(BPI_norm_BMD_select %>%
                       dplyr::rename("bmd" = finalBMD), prob = 0.05)        %>% mutate(method = "p5"),
  get_nth_ranked_bmd(BPI_norm_BMD_select %>%
                       dplyr::rename("bmd" = finalBMD), n = 25)       %>% mutate(method = "rank25"),
  get_LCRD_bmd(BPI_norm_BMD_select %>%
                 dplyr::rename("bmd" = finalBMD))      %>% mutate(method = "LCRD"),
  get_first_mode_bmd(BPI_norm_BMD_select %>%
                       dplyr::rename("bmd" = finalBMD)) %>% mutate(method = "first_mode"),
  get_kneedle_bmd(BPI_norm_BMD_select %>%
                    dplyr::rename("bmd" = finalBMD)) %>% mutate(method = "kneedle")
)



tpod_orig_GU <- bind_rows(
  get_percentile_bmd(GU_norm_BMD_select %>%
                       dplyr::rename("bmd" = finalBMD), prob = 0.05)        %>% mutate(method = "p5"),
  get_nth_ranked_bmd(GU_norm_BMD_select %>%
                       dplyr::rename("bmd" = finalBMD), n = 25)       %>% mutate(method = "rank25"),
  get_LCRD_bmd(GU_norm_BMD_select %>%
                 dplyr::rename("bmd" = finalBMD))      %>% mutate(method = "LCRD"),
  get_first_mode_bmd(GU_norm_BMD_select %>%
                       dplyr::rename("bmd" = finalBMD)) %>% mutate(method = "first_mode"),
  get_kneedle_bmd(GU_norm_BMD_select %>%
                    dplyr::rename("bmd" = finalBMD)) %>% mutate(method = "kneedle")
)



tpod_orig_AU <- bind_rows(
  get_percentile_bmd(AU_norm_BMD_select %>%
                       dplyr::rename("bmd" = finalBMD), prob = 0.05)        %>% mutate(method = "p5"),
  get_nth_ranked_bmd(AU_norm_BMD_select %>%
                       dplyr::rename("bmd" = finalBMD), n = 25)       %>% mutate(method = "rank25"),
  get_LCRD_bmd(AU_norm_BMD_select %>%
                 dplyr::rename("bmd" = finalBMD))      %>% mutate(method = "LCRD"),
  get_first_mode_bmd(AU_norm_BMD_select %>%
                       dplyr::rename("bmd" = finalBMD)) %>% mutate(method = "first_mode"),
  get_kneedle_bmd(AU_norm_BMD_select %>%
                    dplyr::rename("bmd" = finalBMD)) %>% mutate(method = "kneedle")
)





# generate table (NEED TO ADD CI)
tpod_orig_all <- list(
  BMDExpress_log2CPM_noWTT    = tpod_orig_LU,
  BMDExpress_log2CPM_WTT  = tpod_orig_SC,
  DRomics_VST_QT     = tpod_orig_GU,
  DRomics_log2Internal_QT        = tpod_orig_BPI,
  DRomics_CPM_QT = tpod_orig_AU
) %>%
  imap_dfr(~ mutate(.x, analysis_summary = .y)) %>%
  mutate(
    timepoint = factor(timepoint, levels = c("4h","8h","16h","24h","48h","72h")),
    method = factor(method,
                    levels = c("p5","rank25","LCRD","first_mode", "kneedle"),
                    labels = c("5th percentile","25th ranked gene","LCRD","First mode", "Kneedle")
    )
  ) %>%
  rename(tpod_orig = tpod)


data.table::fwrite(tpod_orig_all, file.path(getwd(), "output", "EUT046", "tpod_calculations_alltimepoints.txt"), sep = "\t")



