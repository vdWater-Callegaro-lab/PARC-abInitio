
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

# Kneedle / inflection-point BMD (per timepoint), matching LCRD structure
get_kneedle_bmd <- function(df, S = 1, min_points = 3) {
  stopifnot(is.data.frame(df), all(c("timepoint", "bmd") %in% names(df)))
  
  compute_knee <- function(y_raw, S = 1, min_points = 3) {
    # need enough points to detect an inflection
    y_raw <- y_raw[is.finite(y_raw)]
    if (length(y_raw) < min_points) return(NA_real_)
    
    # Order by BMD; x is rank, y is sorted BMD
    ord <- order(y_raw)
    x   <- seq_along(y_raw)
    y   <- y_raw[ord]
    
    # Smooth the curve (fall back to raw if spline fails)
    ss <- tryCatch(stats::smooth.spline(x, y), error = function(e) NULL)
    if (!is.null(ss)) { x_s <- ss$x; y_s <- ss$y } else { x_s <- x; y_s <- y }
    
    # Normalize x,y to [0,1] (guard zero ranges)
    rx <- range(x_s, finite = TRUE); ry <- range(y_s, finite = TRUE)
    x_n <- if (diff(rx) == 0) rep(0, length(x_s)) else (x_s - rx[1]) / diff(rx)
    y_n <- if (diff(ry) == 0) rep(0, length(y_s)) else (y_s - ry[1]) / diff(ry)
    
    # Difference curve D(x) = y - x
    Dd <- data.frame(x = x_n, y = y_n - x_n, i = seq_along(x_n))
    
    # Find local maxima of D(x)
    if (nrow(Dd) < 3) {
      knee_i <- which.max(Dd$y)
    } else {
      lmx_idx <- which(diff(sign(diff(Dd$y))) == -2) + 1
      if (length(lmx_idx) == 0) {
        knee_i <- which.max(Dd$y)
      } else {
        Dlmx   <- Dd[lmx_idx, , drop = FALSE]
        diff_x <- mean(diff(Dd$x))
        Tlmx   <- Dlmx$y - S * diff_x
        
        # Search regions after each local max for first index below its threshold
        knee_idx <- integer(0)
        n <- nrow(Dd)
        for (k in seq_len(nrow(Dlmx))) {
          start <- Dlmx$i[k]
          end   <- if (k < nrow(Dlmx)) Dlmx$i[k + 1] else n
          idx   <- start:end
          below <- which(Dd$y[idx] < Tlmx[k])
          if (length(below) > 0) knee_idx <- c(knee_idx, idx[below])
        }
        
        knee_i <- if (length(knee_idx) == 0) which.max(Dd$y) else min(knee_idx)
      }
    }
    
    # Map knee index back to the ORIGINAL (unsmoothed) sorted y
    as.numeric(y[knee_i])
  }
  
  df %>%
    dplyr::group_by(timepoint) %>%
    dplyr::summarise(kneedle_bmd = compute_knee(bmd, S = S, min_points = min_points),
                     .groups = "drop")
}




# Summarize BMD metrics per timepoint (simple text table; no knitr/kableExtra)
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
  
  # Compute metrics (using your helpers)
  p5    <- get_percentile_bmd(df, prob = 0.05)        |> dplyr::rename(`P5` = bmd_percentile)
  n25   <- get_nth_ranked_bmd(df, n = 25)             |> dplyr::rename(`Nth25` = nth_bmd)
  lcrd  <- get_LCRD_bmd(df)                           |> dplyr::rename(`LCRD` = lcrd_bmd)
  mode1 <- get_first_mode_bmd(df)                     |> dplyr::rename(`FirstMode` = first_mode_bmd)
  kneedle <- get_kneedle_bmd(df)                      |> dplyr::rename(`Kneedle` = kneedle_bmd)
  
  # Join everything together (keeping all timepoints)
  out <- base |>
    dplyr::left_join(p5,    by = "timepoint") |>
    dplyr::left_join(n25,   by = "timepoint") |>
    dplyr::left_join(lcrd,  by = "timepoint") |>
    dplyr::left_join(mode1, by = "timepoint") |>
    dplyr::left_join(kneedle, by = "timepoint") |>
    dplyr::mutate(timepoint = as.character(timepoint))
  
  # Round numerics for display
  num_cols <- names(out)[sapply(out, is.numeric)]
  out[num_cols] <- lapply(out[num_cols], function(x) round(x, digits))
  
  if (return_data) return(out)
  
  make_kable_tbl <- function(out,
                             digits = 2,
                             caption = NULL) {
    tbl <- out
    if ("timepoint" %in% names(tbl))
      names(tbl)[names(tbl) == "timepoint"] <- "Timepoint"
    
    fmt_num <- function(x)
      ifelse(is.na(x), "NA", sprintf(paste0("%.", digits, "f"), x))
    for (nm in names(tbl)) {
      if (is.numeric(tbl[[nm]])) {
        tbl[[nm]] <- fmt_num(tbl[[nm]])
      } else {
        tbl[[nm]] <- as.character(tbl[[nm]])
        tbl[[nm]][is.na(tbl[[nm]])] <- "NA"
      }
    }
    
    knitr::kable(
      tbl,
      caption = caption,
      align = c("l", rep("r", ncol(tbl) - 1)),
      booktabs = TRUE
    ) |>
      kableExtra::kable_styling(full_width = FALSE)
  }
  
  make_kable_tbl(out, digits = digits, caption = caption)
  
}

