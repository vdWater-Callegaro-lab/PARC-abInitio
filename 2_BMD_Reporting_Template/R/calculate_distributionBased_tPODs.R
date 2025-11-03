
# 5th percentile (or any percentile)
get_percentile_bmd <- function(df, prob = 0.05, na.rm = TRUE) {
  df %>%
    group_by(timepoint) %>%
    summarise(bmd_percentile = quantile(bmd, probs = prob, na.rm = na.rm),
              .groups = "drop")
}

# N-th ranked (e.g., 25th)
get_nth_ranked_bmd <- function(df, n = 25) {
  df %>%
    group_by(timepoint) %>%
    arrange(bmd, .by_group = TRUE) %>%
    mutate(row = row_number()) %>%
    filter(row == n) %>%
    summarise(nth_bmd = first(bmd), .groups = "drop")
}

# Lowest Consistent Ranked Dose (LCRD)
get_LCRD_bmd <- function(df, ratio_threshold = 1.66) {
  compute_lcrd <- function(x) {
    x <- sort(x)
    if (length(x) < 2) return(NA_real_)
    groups <- list(current = x[1])
    largest <- groups$current
    for (i in 2:length(x)) {
      if (x[i] / x[i - 1] <= ratio_threshold) {
        groups$current <- c(groups$current, x[i])
      } else {
        if (length(groups$current) > length(largest)) largest <- groups$current
        groups$current <- x[i]
      }
    }
    if (length(groups$current) > length(largest)) largest <- groups$current
    min(largest)
  }
  
  df %>%
    group_by(timepoint) %>%
    summarise(lcrd_bmd = compute_lcrd(bmd), .groups = "drop")
}

# First mode of the BMD distribution
get_first_mode_bmd <- function(df, bw = "nrd0") {
  first_mode <- function(x) {
    x <- x[is.finite(x)]
    if (length(unique(x)) < 2) return(NA_real_)
    dens <- density(x, bw = bw)
    y <- dens$y
    xgrid <- dens$x
    idx <- which(diff(sign(diff(y))) == -2) + 1  # local maxima
    if (length(idx) == 0) return(NA_real_)
    xgrid[min(idx)]  # first mode
  }
  
  df %>%
    group_by(timepoint) %>%
    summarise(first_mode_bmd = first_mode(bmd), .groups = "drop")
}




# Summarize BMD metrics per timepoint and render a nice Rmd table
bmd_summary_table <- function(
    df,
    digits = 3,
    caption = "BMD summary per timepoint",
    return_data = FALSE
) {
  stopifnot(is.data.frame(df), all(c("timepoint","bmd") %in% names(df)))
  
  # Base timepoint order (natural if order_tp is available)
  tps <- unique(as.character(df$timepoint))
  tps_ord <- tryCatch(levels(order_tp(tps)), error = function(e) tps)
  base <- tibble::tibble(timepoint = factor(tps_ord, levels = tps_ord))
  
  # Compute metrics (using your already-defined helpers)
  p5   <- get_percentile_bmd(df, prob = 0.05)         |> dplyr::rename(`P5` = bmd_percentile)
  n25  <- get_nth_ranked_bmd(df, n = 25)              |> dplyr::rename(`Nth25` = nth_bmd)
  lcrd <- get_LCRD_bmd(df)                            |> dplyr::rename(`LCRD` = lcrd_bmd)
  mode1<- get_first_mode_bmd(df)                      |> dplyr::rename(`FirstMode` = first_mode_bmd)
  
  # Join everything together (keeping all timepoints)
  out <- base |>
    dplyr::left_join(p5,   by = "timepoint") |>
    dplyr::left_join(n25,  by = "timepoint") |>
    dplyr::left_join(lcrd, by = "timepoint") |>
    dplyr::left_join(mode1,by = "timepoint") |>
    dplyr::mutate(timepoint = as.character(timepoint))
  
  # Round numerics for display
  num_cols <- names(out)[sapply(out, is.numeric)]
  out[num_cols] <- lapply(out[num_cols], function(x) round(x, digits))
  
  if (return_data) return(out)
  
  # Pretty table for Rmd
  knitr::kable(
    out,
    caption = caption,
    align = c("l", rep("r", ncol(out) - 1)),
    format = "html"
  ) |>
    kableExtra::kable_styling(
      bootstrap_options = c("striped", "hover", "condensed"),
      full_width = FALSE
    ) |>
    kableExtra::add_header_above(c(" " = 1, "BMD metrics" = ncol(out) - 1))
}

