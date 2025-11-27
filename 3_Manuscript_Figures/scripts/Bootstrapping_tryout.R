

# BOOTSTRAPPING TO IDENTIFY CONFIDENCE INTERVAL IN TPOD CALCULATIONS


## LOAD FUNCTION TO CALUCATE TPODS
source(file.path(getwd(), "R", "calculate_distributionBased_tPODs.R"))


## LOAD DATA
load(file.path(outputDir, "WrangledInput", "WrangledInputData.RData"))



## ORIGINAL TPODS

### CALCULATE ORIGINAL TPODS
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


### COMBINE IN ONE DF
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
    
    # 3) Stack all methods and tag with bootstrap iteration
    bind_rows(res_5, res_25, res_lcrd, res_mode) %>%
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
  LeidenU    = bootstrapped_tpods_LU,
  Sciensano  = bootstrapped_tpods_SC,
  GhentU     = bootstrapped_tpods_GU,
  BPI        = bootstrapped_tpods_BPI,
  AristotleU = bootstrapped_tpods_AU
) %>%
  imap_dfr(~ mutate(.x, partner = .y))   # add partner column from list names

data.table::fwrite(boot_all, file.path(getwd(), "output", "EUT046", "tpod_bootstrapping_results.txt"), sep = "\t")





### GET CONFIDENCE INTERVALS PER METHOD X PARTNER X TIMEPOINT
tpod_ci_all <- boot_all %>%
  group_by(partner, timepoint, method) %>%
  summarise(
    tpod_median = median(tpod, na.rm = TRUE),
    tpod_lower  = quantile(tpod, 0.025, na.rm = TRUE),
    tpod_upper  = quantile(tpod, 0.975, na.rm = TRUE),
    .groups = "drop"
  )


data.table::fwrite(tpod_ci_all, file.path(getwd(), "output", "EUT046", "tpod_bootstrapping_ci_median.txt"), sep = "\t")



## PLOT


### ORDER TIMEPOINT & METHOD FOR PLOT
boot_all <- boot_all %>%
  mutate(
    timepoint = factor(timepoint, levels = c("4h","8h","16h","24h","48h","72h")),
    method = factor(method,
                    levels = c("p5","rank25","LCRD","first_mode"),
                    labels = c("5th percentile","25th ranked gene","LCRD","First mode")
    )
  )



### GENERATE DF FOR PLOTTING
tpod_plot_df <- tpod_ci_all %>%
  left_join(tpod_orig_all,
            by = c("partner","timepoint","method"))



### OPTION 1: FACET BY TPOD METHOD
ggplot(tpod_plot_df,
       aes(x = timepoint, y = tpod_median,
           ymin = tpod_lower, ymax = tpod_upper,
           color = partner)) +
  geom_pointrange(
    position = position_dodge(width = 0.6)
  ) +
  geom_point(
    aes(y = tpod_orig),
    shape = 21,
    fill = "white",
    size = 2.3,
    stroke = 0.6,
    position = position_dodge(width = 0.6)
  ) +
  facet_wrap(~ method, scales = "free_y") +
  labs(
    x = "",
    y = "tPOD",
    color = "",
    title = "tPOD estimates and bootstrap uncertainty"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 12),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_color_manual(values = c("firebrick", "gold1", "gray30", "blue2", "forestgreen"),
                     breaks = c("AristotleU", "GhentU", "BPI", "LeidenU", "Sciensano"))






### OPTION 2: FACET BY TPOD MEASURE AND PARTNER
ggplot(tpod_plot_df,
       aes(x = timepoint, y = tpod_median,
           ymin = tpod_lower, ymax = tpod_upper, color = partner)) +
  geom_pointrange() +
  geom_point(
    aes(y = tpod_orig),
    shape = 21, fill = "white", size = 2.3, stroke = 0.6
  ) +
  facet_grid(method ~ partner, scales = "free_y") +
  labs(
    x = "",
    y = "tPOD",
    title = "Distribution-based tPODs with bootstrap CIs"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_color_manual(values = c("firebrick", "gold1", "gray30", "blue2", "forestgreen"),
                     breaks = c("AristotleU", "GhentU", "BPI", "LeidenU", "Sciensano"))


