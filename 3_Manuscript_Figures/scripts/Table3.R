


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
                       dplyr::rename("bmd" = finalBMD)) %>% mutate(method = "first_mode")
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
                       dplyr::rename("bmd" = finalBMD)) %>% mutate(method = "first_mode")
)



tpod_orig_BPI <- bind_rows(
  get_percentile_bmd(BPI_norm_BMD_select %>%
                       dplyr::rename("bmd" = finalBMD), prob = 0.05)        %>% mutate(method = "p5"),
  get_nth_ranked_bmd(BPI_norm_BMD_select %>%
                       dplyr::rename("bmd" = finalBMD), n = 25)       %>% mutate(method = "rank25"),
  get_LCRD_bmd(BPI_norm_BMD_select %>%
                 dplyr::rename("bmd" = finalBMD))      %>% mutate(method = "LCRD"),
  get_first_mode_bmd(BPI_norm_BMD_select %>%
                       dplyr::rename("bmd" = finalBMD)) %>% mutate(method = "first_mode")
)



tpod_orig_GU <- bind_rows(
  get_percentile_bmd(GU_norm_BMD_select %>%
                       dplyr::rename("bmd" = finalBMD), prob = 0.05)        %>% mutate(method = "p5"),
  get_nth_ranked_bmd(GU_norm_BMD_select %>%
                       dplyr::rename("bmd" = finalBMD), n = 25)       %>% mutate(method = "rank25"),
  get_LCRD_bmd(GU_norm_BMD_select %>%
                 dplyr::rename("bmd" = finalBMD))      %>% mutate(method = "LCRD"),
  get_first_mode_bmd(GU_norm_BMD_select %>%
                       dplyr::rename("bmd" = finalBMD)) %>% mutate(method = "first_mode")
)



tpod_orig_AU <- bind_rows(
  get_percentile_bmd(AU_norm_BMD_select %>%
                       dplyr::rename("bmd" = finalBMD), prob = 0.05)        %>% mutate(method = "p5"),
  get_nth_ranked_bmd(AU_norm_BMD_select %>%
                       dplyr::rename("bmd" = finalBMD), n = 25)       %>% mutate(method = "rank25"),
  get_LCRD_bmd(AU_norm_BMD_select %>%
                 dplyr::rename("bmd" = finalBMD))      %>% mutate(method = "LCRD"),
  get_first_mode_bmd(AU_norm_BMD_select %>%
                       dplyr::rename("bmd" = finalBMD)) %>% mutate(method = "first_mode")
)





# generate table (NEED TO ADD CI)
tpod_orig_all <- list(
  LeidenU    = tpod_orig_LU,
  Sciensano  = tpod_orig_SC,
  GhentU     = tpod_orig_GU,
  BPI        = tpod_orig_BPI,
  AristotleU = tpod_orig_AU
) %>%
  imap_dfr(~ mutate(.x, partner = .y)) %>%
  mutate(
    timepoint = factor(timepoint, levels = c("4h","8h","16h","24h","48h","72h")),
    method = factor(method,
                    levels = c("p5","rank25","LCRD","first_mode"),
                    labels = c("5th percentile","25th ranked gene","LCRD","First mode")
    )
  ) %>%
  rename(tpod_orig = tpod)   # adjust if your col has a different name
