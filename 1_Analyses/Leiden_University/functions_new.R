# functions 

# cpm normalization
cpm_normalization <- function(x) {
  x/sum(x) * 1000000}




# relevance filter
relevance_filter <- function(countdata, metadata) {
  
  # countdata is a dataframe with probe_id as rownames and columns as samples
  # metadata is a dataframe containing sample_id and mean_id columns
  
  countdata[countdata >= 1] <- 1
  countdata[countdata < 1] <- 0
  
  countdata <- melt.data.table(data = data.table(countdata %>%
                                                   rownames_to_column(var = "probe_id")), id.vars = "probe_id", measure.vars = colnames(countdata),
                               variable.name = "SAMPLE_ID", value.name = "count")
  
  countdata <- countdata %>%
    left_join(y = metadata %>%
                select(SAMPLE_ID, TREATMENT_ID), by = "SAMPLE_ID")
  
  # number samples per mean_id group
  output <- countdata %>%
    group_by(probe_id, TREATMENT_ID) %>%
    summarise(count_sum = sum(count)) %>%
    left_join(y = countdata %>%
                distinct(SAMPLE_ID, TREATMENT_ID) %>%
                group_by(TREATMENT_ID) %>%
                mutate(n_TREATMENT_ID = n()) %>%
                distinct(TREATMENT_ID, n_TREATMENT_ID), by = "TREATMENT_ID") %>%
    ungroup() %>%
    mutate(keep = if_else(condition = count_sum >= n_TREATMENT_ID * 0.75, true = TRUE,
                          false = FALSE)) %>%
    group_by(probe_id) %>%
    summarise(keep = sum(keep)) %>%
    mutate(keep = if_else(keep >= 1, true = TRUE, false = FALSE))
  
  return(output)
  
}



# # deseq function
# deseq_function <- function(contrast, countdata, metadata, design = ~ bxSourcePatho, 
#                            shrinkage = FALSE, shrinkage_type = "normal", workers = 1) {
#   
#   # Ensure required packages are available
#   required_packages <- c("DESeq2", "BiocParallel", "dplyr", "tibble")
#   invisible(lapply(required_packages, function(pkg) {
#     if (!requireNamespace(pkg, quietly = TRUE)) {
#       stop(sprintf("Package '%s' is required but not installed. Please install it.", pkg))
#     }
#   }))
#   
#   # Validate the design parameter
#   if(!inherits(design, "formula")) {
#     stop("The 'design' parameter must be a valid formula.")
#   }
#   
#   # Validate shrinkage types
#   valid_shrinkage_types <- c("normal", "apeglm", "ashr")
#   if(shrinkage){
#     if(!(shrinkage_type %in% valid_shrinkage_types)){
#       stop("Invalid `shrinkage_type`. Choose from 'normal', 'apeglm', or 'ashr'.")
#     }
#     # Check if required packages for shrinkage are installed
#     if(shrinkage_type == "apeglm" && !requireNamespace("apeglm", quietly = TRUE)){
#       stop("The 'apeglm' package is required for 'apeglm' shrinkage. Please install it using BiocManager::install('apeglm').")
#     }
#     if(shrinkage_type == "ashr" && !requireNamespace("ashr", quietly = TRUE)){
#       stop("The 'ashr' package is required for 'ashr' shrinkage. Please install it using BiocManager::install('ashr').")
#     }
#   }
#   
#   # Set parallel parameters
#   if (.Platform$OS.type == "windows") {
#     bpparam <- BiocParallel::SnowParam(workers = workers)
#   } else {
#     bpparam <- BiocParallel::MulticoreParam(workers = workers)
#   }
#   
#   # Extract unique pathologies
#   allPathologies <- unique(contrast$main_diagnosis_splitgroups)
#   
#   # Initialize an empty tibble for results
#   deseq_results <- tibble()
#   
#   # Loop over each pathology
#   for (pathology in allPathologies) {
#     message("Processing pathology: ", pathology)
#     
#     # Filter contrasts for the current pathology
#     tmp_contrast <- contrast %>%
#       filter(main_diagnosis_splitgroups == pathology)
#     
#     # Filter metadata for the relevant samples
#     tmp_metadata <- metadata %>%
#       filter(
#         bxSourcePatho %in% c(
#           tmp_contrast$caseCondition,
#           tmp_contrast$controlCondition
#         )
#       )
#     
#     # Select relevant count data
#     tmp_countdata <- countdata %>%
#       select(gene_symbol, all_of(tmp_metadata$sample_id))
#     
#     # Check consistency between metadata and count data
#     if (!all(tmp_metadata$sample_id %in% colnames(tmp_countdata))) {
#       stop("Mismatch between metadata samples and countdata columns for pathology: ", pathology)
#     }
#     
#     # Validate that all design variables are present in metadata
#     design_vars <- all.vars(design)
#     missing_vars <- setdiff(design_vars, colnames(tmp_metadata))
#     if(length(missing_vars) > 0) {
#       stop("The following variables in the design formula are missing from metadata: ", paste(missing_vars, collapse = ", "))
#     }
#     
#     # Create DESeq2 dataset
#     dds <- DESeqDataSetFromMatrix(
#       countData = as.matrix(column_to_rownames(tmp_countdata, var = "gene_symbol")),
#       colData = tmp_metadata,
#       design = design
#     )
#     
#     # Run DESeq2 analysis
#     dds <- DESeq(dds, fitType = "parametric", parallel = TRUE, BPPARAM = bpparam)
#     
#     # Apply each contrast
#     for (j in seq_len(nrow(tmp_contrast))) {
#       contrast_vector <- c(
#         "bxSourcePatho",
#         tmp_contrast$caseCondition[j],
#         tmp_contrast$controlCondition[j]
#       )
#       
#       if(shrinkage){
#         # Apply LFC shrinkage
#         res <- lfcShrink(
#           dds,
#           contrast = contrast_vector,
#           type = shrinkage_type,
#           parallel = TRUE,
#           BPPARAM = bpparam
#         )
#       } else {
#         # Obtain standard results
#         res <- results(
#           object = dds,
#           contrast = contrast_vector,
#           parallel = TRUE,
#           BPPARAM = bpparam
#         )
#       }
#       
#       # Convert results to data frame and add gene symbol
#       res_df <- as.data.frame(res) %>%
#         rownames_to_column("gene_symbol")
#       
#       # If shrinkage is applied, include shrinkage_type
#       if(shrinkage){
#         res_df <- res_df %>%
#           mutate(shrinkage_type = shrinkage_type)
#       }
#       
#       # Add additional columns for context
#       res_df <- res_df %>%
#         mutate(
#           main_diagnosis_splitgroups = pathology,
#           caseCondition = tmp_contrast$caseCondition[j],
#           controlCondition = tmp_contrast$controlCondition[j],
#           shrinkage_applied = shrinkage,
#           shrinkage_type = if_else(shrinkage, shrinkage_type, NA_character_)
#         )
#       
#       # Combine with the results tibble
#       deseq_results <- bind_rows(deseq_results, res_df)
#     }
#   }
#   
#   return(deseq_results)
# }





# Generate DESeq2 dataset and run DESeq
generate_dds <- function(countdata, metadata, contrast, design = ~ MEAN_ID, workers = 1) {
  # Ensure required packages are available
  required_packages <- c("DESeq2", "BiocParallel", "dplyr", "tibble")
  invisible(lapply(required_packages, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required but not installed. Please install it.", pkg))
    }
  }))
  
  # Validate the design parameter
  if(!inherits(design, "formula")) {
    stop("The 'design' parameter must be a valid formula.")
  }
  
  # Set parallel parameters
  if (.Platform$OS.type == "windows") {
    bpparam <- BiocParallel::SnowParam(workers = workers)
  } else {
    bpparam <- BiocParallel::MulticoreParam(workers = workers)
  }
  
  # Extract unique pathologies
  all_experiments = unique(contrast_table$EXPERIMENT)
  
  # Initialize an empty list to store dds objects
  dds_list <- list()
  
  # Loop over each experiment
  for (experiment in all_experiments) {
    message("Processing experiment: ", experiment)
    
    # Filter contrasts for the current experiment
    tmp_contrast <- contrast %>%
      filter(EXPERIMENT == experiment)
    
    # Filter metadata for the relevant samples
    tmp_metadata <- metadata %>%
      filter(
        MEAN_ID %in% c(
          tmp_contrast$MEAN_ID_EXPERIMENT,
          tmp_contrast$MEAN_ID_CONTROL
        )
      )
    
    # Select relevant count data
    tmp_countdata <- countdata %>%
      select(gene_symbol_entrez_id, all_of(tmp_metadata$SAMPLE_ID))
    
    # Check consistency between metadata and count data
    if (!all(tmp_metadata$SAMPLE_ID %in% colnames(tmp_countdata))) {
      stop("Mismatch between metadata samples and countdata columns for experiment: ", experiment)
    }
    
    # Create DESeq2 dataset
    dds <- DESeqDataSetFromMatrix(
      countData = tmp_countdata %>% column_to_rownames("gene_symbol_entrez_id") %>% as.matrix(),
      colData = tmp_metadata,
      design = design
    )
    
    # Run DESeq2 analysis
    dds <- DESeq(dds, fitType = "parametric", parallel = TRUE, BPPARAM = bpparam)
    
    # Store the dds object
    dds_list[[experiment]] <- dds
  }
  
  return(dds_list)
}




# Run DESeq2 analysis with contrasts and shrinkage
run_contrast <- function(dds_list, contrast, shrinkage = FALSE, shrinkage_type = "normal", workers = 1) {
  # Ensure required packages are available
  required_packages <- c("DESeq2", "BiocParallel", "dplyr", "tibble")
  invisible(lapply(required_packages, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required but not installed. Please install it.", pkg))
    }
  }))
  
  # Validate shrinkage types
  valid_shrinkage_types <- c("normal", "apeglm", "ashr")
  if(shrinkage){
    if(!(shrinkage_type %in% valid_shrinkage_types)){
      stop("Invalid `shrinkage_type`. Choose from 'normal', 'apeglm', or 'ashr'.")
    }
    if(shrinkage_type == "apeglm" && !requireNamespace("apeglm", quietly = TRUE)){
      stop("The 'apeglm' package is required for 'apeglm' shrinkage. Please install it using BiocManager::install('apeglm').")
    }
    if(shrinkage_type == "ashr" && !requireNamespace("ashr", quietly = TRUE)){
      stop("The 'ashr' package is required for 'ashr' shrinkage. Please install it using BiocManager::install('ashr').")
    }
  }
  
  # Set parallel parameters
  if (.Platform$OS.type == "windows") {
    bpparam <- BiocParallel::SnowParam(workers = workers)
  } else {
    bpparam <- BiocParallel::MulticoreParam(workers = workers)
  }
  
  # Initialize an empty tibble for results
  deseq_results <- tibble()
  
  # Extract unique experiments
  all_experiments <- names(dds_list)
  
  # Loop over each experiment
  for (experiment in all_experiments) {
    message("Processing experiment: ", experiment)
    
    dds <- dds_list[[experiment]]
    
    # Filter contrasts for the current experiment
    tmp_contrast <- contrast %>%
      filter(EXPERIMENT == experiment)
    
    # Apply each contrast
    for (j in seq_len(nrow(tmp_contrast))) {
      contrast_vector <- c(
        "MEAN_ID",
        tmp_contrast$MEAN_ID_EXPERIMENT[j],
        tmp_contrast$MEAN_ID_CONTROL[j]
      )
      
      if(shrinkage){
        # Apply LFC shrinkage
        res <- lfcShrink(
          dds,
          contrast = contrast_vector,
          type = shrinkage_type,
          parallel = TRUE,
          BPPARAM = bpparam
        )
      } else {
        # Obtain standard results
        res <- results(
          object = dds,
          contrast = contrast_vector,
          parallel = TRUE,
          BPPARAM = bpparam
        )
      }
      
      # Convert results to data frame and add gene symbol
      res_df <- as.data.frame(res) %>%
        rownames_to_column("gene_symbol_entrez_id")
      
      
      # Add additional columns for context
      res_df <- res_df %>%
        mutate(
          experiment = experiment,
          MEAN_ID_experiment = tmp_contrast$MEAN_ID_EXPERIMENT[j],
          MEAN_ID_control = tmp_contrast$MEAN_ID_CONTROL[j],
          shrinkage_applied = shrinkage,
          shrinkage_type = if_else(shrinkage, shrinkage_type, NA_character_)
        )
      
      # Combine with the results tibble
      deseq_results <- bind_rows(deseq_results, res_df)
    }
  }
  
  return(deseq_results)
}