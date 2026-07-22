
studyNR <- "EUT046"


library(tidyverse)
library(data.table)

metadata <- fread("data/EUT046/EUT046_metadata_final.csv") 

timepoints_to_permute <- metadata %>% pull(TIMEPOINT) %>% unique()

set.seed(123)


make_permuted_metadata <- function(meta, perm_id, seed) {
  set.seed(seed)
  
  meta_perm <- meta %>%
    mutate(
      original_concentration = CONCENTRATION,
      permuted_concentration = CONCENTRATION
    )
  
  for (tp in timepoints_to_permute) {
    idx <- meta_perm$TIMEPOINT == tp
    
    meta_perm$permuted_concentration[idx] <- sample(
      meta_perm$CONCENTRATION[idx],
      size = sum(idx),
      replace = FALSE
    )
  }
  
  meta_perm <- meta_perm %>%
    mutate(
      permutation_id = perm_id,
      permuted_treatment = if_else(
        permuted_concentration == 0,
        "control",
        "cisplatin"
      )
    ) %>%
    group_by(TIMEPOINT, permuted_concentration) %>%
    mutate(permuted_replicate = row_number()) %>%
    ungroup()
  
  return(meta_perm)
}

check_permutation <- function(meta_original, meta_perm) {
  original_counts <- meta_original %>%
    filter(TIMEPOINT %in% timepoints_to_permute) %>%
    count(TIMEPOINT, CONCENTRATION, name = "n_original")
  
  permuted_counts <- meta_perm %>%
    filter(TIMEPOINT %in% timepoints_to_permute) %>%
    count(TIMEPOINT, permuted_concentration, name = "n_permuted") %>%
    rename(CONCENTRATION = permuted_concentration)
  
  check <- left_join(
    original_counts,
    permuted_counts,
    by = c("TIMEPOINT", "CONCENTRATION")
  )
  
  all(check$n_original == check$n_permuted)
}

dir.create(file.path("output", studyNR, "permuted_metadata"), showWarnings = FALSE)

n_perm <- 100
seeds <- 1000 + seq_len(n_perm)

manifest <- tibble(
  permutation_id = sprintf("perm_%03d", seq_len(n_perm)),
  seed = seeds,
  file = file.path(
    "output",
    studyNR,
    "permuted_metadata",
    paste0("metadata_perm_", sprintf("%03d", seq_len(n_perm)), ".csv")
  )
)

for (i in seq_len(n_perm)) {
  perm_id <- manifest$permutation_id[i]
  seed <- manifest$seed[i]
  
  meta_perm <- make_permuted_metadata(metadata, perm_id, seed)
  
  stopifnot(check_permutation(metadata, meta_perm))
  
  write_csv(meta_perm, manifest$file[i])
}

write_csv(manifest, file.path("output", studyNR, "permutation_manifest.csv"))