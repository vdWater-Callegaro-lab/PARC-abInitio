library(DESeq2)
setwd("C:/Users/EmDe5048/Desktop")
#Data importation (RODAF-filtered and summed counts + metadata)
EUT046 <- read.csv("../Desktop/PARC_project/Sys_Tox_project/abInitio/datasets/EUT046_RODAFpreprocessed_sumcounts.csv", row.names = 'gene_symbol')
EUT046_metadata <- read.csv("../Desktop/PARC_project/Sys_Tox_project/abInitio/metadata/EUT046_metadata_desing.csv", row.names = 'SAMPLE_ID')

#Checking if colnames and rownames are identical and in the same order
all(rownames(EUT046_metadata) %in% colnames(EUT046))
all(rownames(EUT046_metadata) == colnames(EUT046))
EUT046 <- EUT046[, rownames(EUT046_metadata)]
all(rownames(EUT046_metadata) == colnames(EUT046))

#Factorization of variables included in design
EUT046_metadata$REPLICATE <- factor(EUT046_metadata$REPLICATE)
EUT046_metadata$DOSE_LEVEL <-  factor(EUT046_metadata$DOSE_LEVEL)
EUT046_metadata$TIME <-  factor(EUT046_metadata$TIME)
EUT046_metadata$TREATMENT <- factor(EUT046_metadata$TREATMENT)
EUT046_metadata$treat_conc_time <-  factor(EUT046_metadata$treat_conc_time)
levels(EUT046_metadata$REPLICATE)
#Design creation
design <- ~ REPLICATE + treat_conc_time

#DDS object building 
dds <- DESeqDataSetFromMatrix(EUT046, EUT046_metadata, design = design)

#Running DESeq2
dds <- DESeq(dds)
EUT046_normalized <- counts(dds, normalized = TRUE)
#write.csv(EUT046_normalized, file = "../Desktop/PARC_project/Sys_Tox_project/abInitio/datasets/EUT046_RODAFpreprocessed_normalizedcounts.csv")

#Isolate treatment and control names to loop on and create contrasts results
sample_names <- levels(EUT046_metadata$treat_conc_time)
treated_names <- grep(pattern = 'CSP', sample_names, value=TRUE)
control_names <- grep(pattern = 'DMEM', sample_names, value=TRUE)

for (e in treated_names) {
  e_splitted <- strsplit(e, '_')
  e_time <- paste(e_splitted[[1]][2], '_', sep = '')
  c <-  paste('DMEM_0_', tail(e_splitted[[1]], n = 1), sep = '')
  print(e)
  print(c)
  res <- results(dds, contrast = c('treat_conc_time', e, c))
  print('done')
  write.csv(res, file = paste("../Desktop/PARC_project/Sys_Tox_project/abInitio/DESeq2_res/EUT046/", e, "_", c, ".csv", sep = ""))
}