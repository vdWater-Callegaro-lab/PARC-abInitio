library(edgeR)
library(ggplot2)
library(gplots)
library(vegan)
library(lme4)
library(emmeans)
library(ggpubr)
library(DESeq2)
library(ape)
library(ballgown)
library(dplyr)
library(tidyr)
library(tibble)
library(magrittr)
library(openxlsx)
library(EnhancedVolcano)

setwd("C:/Users/toxlab/Lab/Benakeio/Kiki/PARC/PARC_5.3.1/data/Ab_Initio/RPTEC-TERT1_CIS")

rawcount.d <- read.csv("EUT046_counts_per_gene_per_sample_raw_final.csv", row.names = "Column1")
rawcount.d$Name <- as.factor(sapply(strsplit(rownames(rawcount.d), split = "_"), "[", 1)) # add Name column
rawcount.dS <- rawcount.d
rawcount.dS <- aggregate(. ~ Name, rawcount.dS, sum)                                      # sum up the genes with the same name
rownames(rawcount.dS) <- rawcount.dS$Name
rawcount.dS <- rawcount.dS[ ,-1]
rawcount.dS$rowSums <- rowSums(rawcount.dS)                                              # create rowSums columm,
rawcount.dS <- rawcount.dS[order(rawcount.dS$rowSums, decreasing = TRUE),]               # and order by Sums column (descending in counts) and...
# rawcount.dS <- rawcount.dS[order(rownames(rawcount.dS), decreasing = FALSE),]            # .... or, order by rowname
rawcount.dS <- rawcount.dS[ ,-dim(rawcount.dS)[2]]                                       # and remove last (rowSums) column

meta.EUT046 <- read.csv("EUT046_metadata_final.csv")
meta.EUT046$SAMPLE <- gsub(" ", "", meta.EUT046$SAMPLE)
meta.EUT046$TREATMENT_ID <- gsub(" ", "", meta.EUT046$TREATMENT_ID)
meta.EUT046[sapply(meta.EUT046, is.character)] <- lapply(meta.EUT046[sapply(meta.EUT046, is.character)], as.factor) # all chr to factors
meta.EUT046[sapply(meta.EUT046, is.numeric)] <- lapply(meta.EUT046[sapply(meta.EUT046, is.numeric)], as.factor) # all num to factors
meta.EUT046 <- meta.EUT046[match(colnames(rawcount.dS), meta.EUT046$SAMPLE_ID),]                                 # reorder by count file column (sample) names

all(meta.EUT046$SAMPLE_ID %in% colnames(rawcount.dS))


dds <- DESeqDataSetFromMatrix(countData = rawcount.dS,
                              colData = meta.EUT046,
                              design= ~ REPLICATE + TREATMENT_ID,
                              tidy = FALSE)

smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep, ]

dds <- DESeq(dds)
resultsNames(dds)



# Contrasts
# see v1 script
#

### while loop for res.* objects to export result tables
i = 1
temp <- objects(pattern = "res\\.CS")

while(i <= length(temp)){
  s <- get(temp[i])
  cat(temp[i])
  cat("\n")
  
  degs <- s[which(s$padj < 0.01),]
  write.table(s, file = paste0(temp[i], "_DEGs.", "txt"))
  
  i = i + 1
}



##### while loop for res.* objects to count the number of significant DEGs and other preliminary operations_________________________________________ START
i = 1
temp <- objects(pattern = "res\\.CS")
sig.lengths <- c()
sig.names <- list()

while(i <= length(temp)){
  s <- get(temp[i])
  cat(temp[i])
  cat("\n")
  sl <- length(which(s$padj<0.05))
  sn <- rownames(s[which(s$padj<0.05),])
  
  sig.lengths[i] <- sl
  sig.names[[i]] <- sn
  i = i + 1
}
names(sig.names) <- temp

names(sig.lengths) <- temp
names(sig.lengths) <- gsub("_8_", "_08_", names(sig.lengths)) # rename 8 to 08 for appropriate ordering
names(sig.lengths) <- gsub("_4_", "_04_", names(sig.lengths)) # rename 4 to 04 for appropriate ordering
names(sig.lengths) <- gsub("_1v0", "_01v0", names(sig.lengths))
names(sig.lengths) <- gsub("_2.5v0", "_02.5v0", names(sig.lengths))
names(sig.lengths) <- gsub("_5v0", "_05v0", names(sig.lengths))
sig.lengths <- sig.lengths[order(names(sig.lengths))] # order sig.lengths

barplot(sig.lengths, col=rainbow(length(sig.lengths)), las=2)

df.sig.lengths <- data.frame(sig.lengths)
df.sig.lengths$timepoint <- sapply(strsplit(rownames(df.sig.lengths), "_"), "[", 2 )
df.sig.lengths$comparison <- sapply(strsplit(rownames(df.sig.lengths), "_"), "[", 3 )
ggplot(df.sig.lengths, aes(x = timepoint, y = sig.lengths, fill = comparison)) + geom_bar(stat = "identity") + ggtitle("EUT046")

# ggsave("NumberOfSignificant.png", width = 4.5, height = 3, dpi = 300)


# Retrieve the normalized count matrix
ddsNC <- counts(dds, normalized = TRUE)

# Log2 transform the normalized counts
ddsNCL2 <- log2(ddsNC + 1)

# Apply VST transformation
vsd_dds <- vst(dds, blind = FALSE)
assay(vsd_dds)

# Compute counts per million (CPM)
cpm_dds <- cpm(counts(dds, normalized = TRUE))

# Apply rlog transformation
rld_dds <- rlog(dds, blind = FALSE)
assay(rld_dds)
se <- SummarizedExperiment(log2(counts(dds, normalized=TRUE) + 1), colData=colData(dds))
plotPCA.san <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC3 = pca$x[, 3], PC4 = pca$x[, 4], group = group, 
                  intgroup.df, name = colData(se)[,1])
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC3", y = "PC4", color = "group", label = "name")) + geom_point(size = 3) + xlab(paste0("PC3: ", round(percentVar[3] * 100), "% variance")) + ylab(paste0("PC4: ", round(percentVar[4] * 100), "% variance")) + coord_fixed() + geom_text_repel(size=3) 
  
}
plotPCA(DESeqTransform(se), "CONCENTRATION")
plotPCA(DESeqTransform(se), "REPLICATE")
plotPCA(DESeqTransform(se), "TREATMENT_ID")
#  ggsave("PCA_EUT046_DOSE_LEVEL.png", width = 9, height = 6, dpi = 300)


# Filter to have at least one condition with 75% of the samples above 1 CPM
# Check: https://github.com/R-ODAF/Main/blob/main/scripts/R_ODAF_DEGs.R
Groups <- table(meta.EUT046[,"TREATMENT_ID"])

Filter <- matrix(data = NA, ncol = 3, nrow = nrow(ddsNC))
rownames(Filter) <- rownames(ddsNC)
colnames(Filter) <- c("Low","quantile","spike")

for (gene in 1:nrow(dds)){
  CountsPass <- NULL
  for (group in 1:length(Groups)){
    sampleCols <- grep(dimnames(Groups)[[1]][group], meta.EUT046[,"TREATMENT_ID"])
    Check <- sum(cpm_dds[gene, sampleCols] >= 0.5) >= 0.75 * Groups[group]
    CountsPass <- c(CountsPass, Check)
  }
  if (sum(CountsPass) > 0) {Filter[gene,1] <- 1}	else {Filter[gene,1] <- 0}
}

ddsNCF <- ddsNC[Filter[,1] == 1,]
ddsF <- dds[Filter[,1] == 1,]

##### while loop for res.* objects to count the number of significant DEGs, Filtering and other preliminary operations_________________________________________END


#####################
# Feature Selection #
#####################
library(diptest)

#### DIP test to the whole data set
# Initialize a vector to store dip test p-values
dip_p_values <- numeric(nrow(ddsNCF))

# Apply the Dip Test to each gene
for (i in 1:nrow(ddsNCF)) {
  gene_expr <- ddsNCF[i, ]
  dip_p_values[i] <- dip.test(gene_expr)$p.value
}

# Adjust p-values for multiple testing (e.g., using Benjamini-Hochberg method)
dip_p_values_adjusted <- p.adjust(dip_p_values, method = "BH")

# Identify significant genes
signif_dip <- rownames(ddsNCF)[dip_p_values_adjusted < 0.05]
print(signif_dip)
# intersect(signif_dip, Reduce(intersect, sig.names[grep("_24_", names(sig.names))]))    # with all common genes at 24h
# intersect(signif_dip, Reduce(intersect, sig.names[grep("_72_", names(sig.names))]))    # with all common genes at 72h

plot(density(ddsNCF[which(rownames(ddsNCF) == "CHEK1"),])); rug(ddsNCF[which(rownames(ddsNCF) == "CHEK1"),])
plot(density(ddsNCF[which(rownames(ddsNCF) == "TIPARP"),])); rug(ddsNCF[which(rownames(ddsNCF) == "TIPARP"),])



#### DIP test by condition
# by TIMEPOINT
conditions <- unique(meta.EUT046$TIMEPOINT)
dip_p_values_byTIMEPOINT <- list()
for (TIMEPOINT in conditions) {
  condition_samples <- meta.EUT046$TIMEPOINT == TIMEPOINT
  dip_p_val <- numeric(nrow(ddsNCF))
  
  for (i in 1:nrow(ddsNCF)) {
    gene_expression <- ddsNCF[i, condition_samples]
    dip_p_val[i] <- dip.test(gene_expression)$p.value
  }
  
  dip_p_values_byTIMEPOINT[[TIMEPOINT]] <- dip_p_val
}

# Adjust p-values for multiple testing within each condition
dip_p_values_byTIMEPOINT_adjusted <- lapply(dip_p_values_byTIMEPOINT, function(p_values) p.adjust(p_values, method = "BH"))

# Identify significant genes for each condition
signif_dip_byTIMEPOINT <- lapply(dip_p_values_byTIMEPOINT_adjusted, function(p_values) rownames(ddsNCF)[p_values < 0.05])

# Print significant genes for each condition
print(signif_dip_byTIMEPOINT)

plot(density(ddsNC[which(rownames(ddsNCF) == "HIST1H3H"),])); rug(ddsNCF[which(rownames(ddsNCF) == "HIST1H3H"),])
plot(density(ddsNC[which(rownames(ddsNCF) == "CDK1"),])); rug(ddsNCF[which(rownames(ddsNCF) == "CDK1"),])


# by DOSE LEVEL
conditions <- unique(meta.EUT046$DOSE_LEVEL)
dip_p_values_byDOSE <- list()
for (DOSE_LEVEL in conditions) {
  condition_samples <- meta.EUT046$DOSE_LEVEL == DOSE_LEVEL
  dip_p_val <- numeric(nrow(ddsNCF))
  
  for (i in 1:nrow(ddsNCF)) {
    gene_expression <- ddsNCF[i, condition_samples]
    dip_p_val[i] <- dip.test(gene_expression)$p.value
  }
  
  dip_p_values_byDOSE[[DOSE_LEVEL]] <- dip_p_val
}

# Adjust p-values for multiple testing within each condition
dip_p_values_byDOSE_adjusted <- lapply(dip_p_values_byDOSE, function(p_values) p.adjust(p_values, method = "BH"))

# Identify significant genes for each condition
signif_dip_byDOSE <- lapply(dip_p_values_byDOSE_adjusted, function(p_values) rownames(ddsNCF)[p_values < 0.05])

# Print significant genes for each condition
print(signif_dip_byDOSE)

plot(density(ddsNCF[which(rownames(ddsNCF) == "HIST1H3F"),])); rug(ddsNCF[which(rownames(ddsNCF) == "HIST1H3F"),])



# by DOSE LEVEL at 72h
ddsNCF_72h <- ddsNCF[ ,grep("_72h_", colnames(ddsNCF))]
meta.EUT046_72h <- meta.EUT046[meta.EUT046$TIMEPOINT == "72hr", ]

conditions <- unique(meta.EUT046_72h$DOSE_LEVEL)
dip_p_values_byDOSE <- list()
for (DOSE_LEVEL in conditions) {
  condition_samples <- meta.EUT046_72h$DOSE_LEVEL == DOSE_LEVEL
  dip_p_val <- numeric(nrow(ddsNCF_72h))
  
  for (i in 1:nrow(ddsNCF_72h)) {
    gene_expression <- ddsNCF_72h[i, condition_samples]
    dip_p_val[i] <- dip.test(gene_expression)$p.value
  }
  
  dip_p_values_byDOSE[[DOSE_LEVEL]] <- dip_p_val
}

# Adjust p-values for multiple testing within each condition
dip_p_values_byDOSE_adjusted <- lapply(dip_p_values_byDOSE, function(p_values) p.adjust(p_values, method = "BH"))

# Identify significant genes for each condition
signif_dip_byDOSE <- lapply(dip_p_values_byDOSE_adjusted, function(p_values) rownames(ddsNCF_72h)[p_values < 0.05])

# Print significant genes for each condition
print(signif_dip_byDOSE)

plot(density(ddsNCF_72h[which(rownames(ddsNCF_72h) == "WBSCR22"),])); rug(ddsNCF_72h[which(rownames(ddsNCF_72h) == "WBSCR22"),])






####################
### BMD Analysis ###
####################
library(DRomics)

DRdose <- as.numeric(gsub("uM", "", meta.EUT046$CONCENTRATION))
df_DRinit <- counts(dds, normalized = FALSE)
timeLevs <- as.character(unique(lapply(strsplit(colnames(df_DRinit), "_"), "[", 2)))[-1]

out.DRomic <- list()
for(i in 1:length(timeLevs)){
  
  print(timeLevs[i])
  df_temp <- df_DRinit[, grep(timeLevs[i], colnames(df_DRinit))]
  dr.dds <- formatdata4DRomics(df_temp, dose = DRdose[grep(timeLevs[i], colnames(df_DRinit))], samplenames = colnames(df_temp))
  df_forDR <- RNAseqdata(dr.dds)
  # selection of significantly responding items
  item.DRomic <- itemselect(df_forDR, select.method = "quadratic", FDR = 0.01)
  # choice of the best fit for each curve
  dr.fit <- drcfit(item.DRomic, progressbar = TRUE)
  #calculation of bmd
  dr.calcBMD <- bmdcalc(dr.fit, z = 1, x = 10)
  # Calculation of confidence intervals on the BMDs by bootstrap
  dr.b <- bmdboot(dr.calcBMD, niter = 1000, progressbar = TRUE, parallel = "snow", ncpus = 12)
  # Filtering of BMD calculation
  dr.filtBMD <- bmdfilter(dr.b$res, BMDtype = "zSD", BMDfilter = "finiteCI")
  
  out.DRomic[[i]] <- list(fit = dr.fit,
                          calcBMD = dr.calcBMD,
                          bootCI = dr.b,
                          filtBMD = dr.filtBMD)
  
}
names(out.DRomic) <- timeLevs


# create stat accounts and outputs (Tables & Plots)
library(ModEstM)

StatRep <- list()
for (i in 1:length(out.DRomic)){
  
  # Calculate the 5th percentile and the 25th lowest ranked gene
  p5 <- quantile(out.DRomic[[i]]$filtBMD$BMD.zSD, 0.05)
  m1 <- ModEstM(out.DRomic[[i]]$filtBMD$BMD.zSD)
  lrg25 <- sort(out.DRomic[[i]]$filtBMD$BMD.zSD)[25]
  
  print(p5)
  print(m1)
  print(lrg25)
  
  StatRep[[i]] <- list(perc5th = p5,
                       mode1st = m1,
                       rank25th = lrg25)
  
  # Export fit results
  # write.table(out.DRomic[[i]]$fit$fitres, file = paste0("Out.EUT046.dr.fit.", timeLevs[i], ".txt"))
  # # Export plot of fitted curves
  # jpeg(paste0("Out.EUT046.dr.fit.", timeLevs[i], ".jpg"), width = 9.6, height = 5.7, units = "in",res = 300, pointsize = 12, quality = 75)
  # plot(out.DRomic[[i]]$fit, dose_log_transfo = FALSE)
  # dev.off()
  # 
  # # Export bmd results
  # write.table(out.DRomic[[i]]$calcBMD$res, file = paste0("Out.EUT046.dr.fit.", timeLevs[i], ".txt"), row.names = FALSE)
  # # Export plot of the BMD distribution
  # jpeg(paste0("Out.EUT046.dr.calcBMD.", timeLevs[i], ".jpg"), width = 9.6, height = 5.7, units = "in",res = 300, pointsize = 12, quality = 75)
  # plot(out.DRomic[[i]]$calcBMD, BMDtype = "zSD", plottype = "ecdf") + theme_bw()
  # dev.off()
  # 
  # # Export unfiltered and filtered bootstrap results
  # write.table(out.DRomic[[i]]$bootCI$res, file = paste0("Out.EUT046.dr.bootCI.", timeLevs[i], ".txt"), row.names = FALSE)
  # write.table(out.DRomic[[i]]$filtBMD, file = paste0("Out.EUT046.dr.filtBMD.", timeLevs[i], ".txt"), row.names = FALSE)
  # 
  # # Export BMD plots
  # jpeg(paste0("Out.EUT046.FiltBMD.", timeLevs[i], ".jpg"), width = 9.6, height = 5.7, units = "in",res = 300, pointsize = 12, quality = 75)
  # bmdplot(out.DRomic[[i]]$filtBMD, BMDtype = "zSD", point.size = 2, point.alpha = 0.4, add.CI = TRUE, line.size = 0.4) + theme_bw()
  # dev.off()
  # # Export the fitted curves with BMD and CI
  # jpeg(paste0("Out.EUT046.FiltBMDbootCI.", timeLevs[i], ".jpg"), width = 9.6, height = 5.7, units = "in",res = 300, pointsize = 12, quality = 75)
  # plot(out.DRomic[[i]]$fit, BMDoutput = out.DRomic[[i]]$bootCI, dose_log_transfo = FALSE)
  # dev.off()
  
}
names(StatRep) <- timeLevs

# Create Accumulation Plot
df_All.Out.DRomic <- bind_rows(list(out.DRomic[[1]]$filtBMD,
                                    out.DRomic[[2]]$filtBMD,
                                    out.DRomic[[3]]$filtBMD,
                                    out.DRomic[[4]]$filtBMD,
                                    out.DRomic[[5]]$filtBMD,
                                    out.DRomic[[6]]$filtBMD
),
.id = "timepoint")

df_All.Out.DRomic$timepoint[which(df_All.Out.DRomic$timepoint == "1")] <- "16h"
df_All.Out.DRomic$timepoint[which(df_All.Out.DRomic$timepoint == "2")] <- "24h"
df_All.Out.DRomic$timepoint[which(df_All.Out.DRomic$timepoint == "3")] <- "48h"
df_All.Out.DRomic$timepoint[which(df_All.Out.DRomic$timepoint == "4")] <- "4h"
df_All.Out.DRomic$timepoint[which(df_All.Out.DRomic$timepoint == "5")] <- "72h"
df_All.Out.DRomic$timepoint[which(df_All.Out.DRomic$timepoint == "6")] <- "8h"



df_All.Out.DRomic$timepoint <- as.factor(df_All.Out.DRomic$timepoint)

BMDaccumulationPlot <- function(BMDoutput, level) {
  hover <- if (level == "gene") "gene: " else "pathway: "
  BMDoutput %>%
    group_by(timepoint) %>%
    arrange(BMD.zSD) %>%
    mutate(cumulative_BMD = cumsum(BMD.zSD)) %>%
    ggplot(aes(BMD.zSD, y = cumulative_BMD, color = timepoint, shape = timepoint)) +
    geom_point(stat = "identity") +
    scale_x_log10() +
    scale_color_discrete(name = "Timepoint", labels = timeLevs) +
    labs(title = "BMD accumulation plot",
         y = "Cumulative BMD", 
         x = "BMD")
}

# jpeg("Accum.Plot.Out.EUT046filtBMD.jpg", width = 9.6, height = 5.7, units = "in",res = 300, pointsize = 12, quality = 75)
BMDaccumulationPlot(BMDoutput = df_All.Out.DRomic, level = "id")
# dev.off()




########################
### Pathway Analysis ###
########################

# https://rgd.mcw.edu/wg/home
# https://yulab-smu.top/biomedical-knowledge-mining-book/dose-enrichment.html
# https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html
library(ReactomePA)
library(pathview)
library(clusterProfiler)
library(DOSE)
library(AnnotationHub)
library(msigdbr)

## Prepare MeSH hub, data frames. match manif(ests) with results tables, etc_______________________________________________________________________________________________ START
ah <- AnnotationHub(localHub=TRUE)
hsa <- query(ah, c("MeSHDb", "Homo sapiens"))
file_hsa <- hsa[[1]]
db <- MeSHDbi::MeSHDb(file_hsa)
hall_set <- msigdbr(species = "Homo sapiens", category = "H")

# "C:\\Users\\toxlab\\AppData\\Local/R/cache/R/AnnotationHub/3b9057687141_98354" 

manif1 <- read.table("C:/Users/toxlab/Lab/Benakeio/Kiki/PARC/PARC_5.3.1/data/Ab_Initio/Manifests/Temposeq_manifest_Human_WT_2.0_realignment_96percent_2023-05-15.txt", 
                     sep = "\t", header = TRUE)
manif1 <- manif1[, c(1,2,3,9,10)]
manif2 <- read.table("C:/Users/toxlab/Lab/Benakeio/Kiki/PARC/PARC_5.3.1/data/Ab_Initio/Manifests/Temposeq_manifest_Human_WT_1.2_realignment_96percent_2023-05-15.txt", 
                     sep = "\t", header = TRUE)
manif2 <- manif2[, c(1,2,3,9,8)]
colnames(manif2) <- colnames(manif1)

manif <- merge(manif1, manif2, all = TRUE)




### For DRomics: ###

# Create list of data frames    *** Careful: Timepoints are ordered in DRomix.framesL (unlike df_All.Out.DRomic)
DRomix.framesL <- list(T4h = out.DRomic[[4]]$filtBMD,
                       T8h = out.DRomic[[6]]$filtBMD,
                       T16h = out.DRomic[[1]]$filtBMD,
                       T24h = out.DRomic[[2]]$filtBMD,
                       T48h = out.DRomic[[3]]$filtBMD,
                       T72h = out.DRomic[[5]]$filtBMD)


# while loop for Gene Set Enrichment Analysis on the DRomix.framesL_________________________________________ START
i = 1
temp <- length(DRomix.framesL)
AllPathDR <- list()

while(i <= temp){
  
  cat(paste("Processing", i, "of", length(DRomix.framesL), "\n", sep = " "))
  s <- DRomix.framesL[[i]]
  df <- data.frame(s)
  gene_symbol <- manif[match(df$id, manif$Gene_Symbol_TS),"gene_symbol"]
  entrez_id <- manif[match(df$id, manif$Gene_Symbol_TS),"Entrez.ID"]
  df <- cbind(df, gene_symbol, entrez_id)
  
  sdrx <- df$entrez_id             # for over-representation analysis (enrichXXX)
  
  
  REACT.enr <- enrichPathway(sdrx)
  KEGG.enr <- enrichKEGG(sdrx)
  HALL.enr <- enricher(sdrx, TERM2GENE = hall_set[,c("gs_name","entrez_gene")], pvalueCutoff = 0.05, pAdjustMethod = "fdr")
  WP.enr <- enrichWP(sdrx, organism = "Homo sapiens")
  
  AllPathDR[[i]] <- list(REA = REACT.enr,
                         KEGG = KEGG.enr,
                         HALL = HALL.enr,
                         WP = WP.enr)
  
  i = i + 1
}
names(AllPathDR) <- names(DRomix.framesL)

## Compile result tables
# REACTOME
df_Path.Out.React <- bind_rows(list(cbind(AllPathDR[[1]]$REA@result, "Timepoint" = rep("4h", dim(AllPathDR[[1]]$REA@result)[1])),
                                    cbind(AllPathDR[[2]]$REA@result, "Timepoint" = rep("8h", dim(AllPathDR[[2]]$REA@result)[1])),
                                    cbind(AllPathDR[[3]]$REA@result, "Timepoint" = rep("16h", dim(AllPathDR[[3]]$REA@result)[1])),
                                    cbind(AllPathDR[[4]]$REA@result, "Timepoint" = rep("24h", dim(AllPathDR[[4]]$REA@result)[1])),
                                    cbind(AllPathDR[[5]]$REA@result, "Timepoint" = rep("48h", dim(AllPathDR[[5]]$REA@result)[1])),
                                    cbind(AllPathDR[[6]]$REA@result, "Timepoint" = rep("72h", dim(AllPathDR[[6]]$REA@result)[1]))
))
df_Path.Out.React$Timepoint <- as.factor(df_Path.Out.React$Timepoint)

# KEGG
df_Path.Out.KEGG <- bind_rows(list(cbind(AllPathDR[[1]]$KEGG@result, "Timepoint" = rep("4h", dim(AllPathDR[[1]]$KEGG@result)[1])),
                                   cbind(AllPathDR[[2]]$KEGG@result, "Timepoint" = rep("8h", dim(AllPathDR[[2]]$KEGG@result)[1])),
                                   cbind(AllPathDR[[3]]$KEGG@result, "Timepoint" = rep("16h", dim(AllPathDR[[3]]$KEGG@result)[1])),
                                   cbind(AllPathDR[[4]]$KEGG@result, "Timepoint" = rep("24h", dim(AllPathDR[[4]]$KEGG@result)[1])),
                                   cbind(AllPathDR[[5]]$KEGG@result, "Timepoint" = rep("48h", dim(AllPathDR[[5]]$KEGG@result)[1])),
                                   cbind(AllPathDR[[6]]$KEGG@result, "Timepoint" = rep("72h", dim(AllPathDR[[6]]$KEGG@result)[1]))
))
df_Path.Out.KEGG$Timepoint <- as.factor(df_Path.Out.KEGG$Timepoint)

# HALL
df_Path.Out.HALL <- bind_rows(list(cbind(AllPathDR[[1]]$HALL@result, "Timepoint" = rep("4h", dim(AllPathDR[[1]]$HALL@result)[1])),
                                   cbind(AllPathDR[[2]]$HALL@result, "Timepoint" = rep("8h", dim(AllPathDR[[2]]$HALL@result)[1])),
                                   cbind(AllPathDR[[3]]$HALL@result, "Timepoint" = rep("16h", dim(AllPathDR[[3]]$HALL@result)[1])),
                                   cbind(AllPathDR[[4]]$HALL@result, "Timepoint" = rep("24h", dim(AllPathDR[[4]]$HALL@result)[1])),
                                   cbind(AllPathDR[[5]]$HALL@result, "Timepoint" = rep("48h", dim(AllPathDR[[5]]$HALL@result)[1])),
                                   cbind(AllPathDR[[6]]$HALL@result, "Timepoint" = rep("72h", dim(AllPathDR[[6]]$HALL@result)[1]))
))
df_Path.Out.HALL$Timepoint <- as.factor(df_Path.Out.HALL$Timepoint)

# WP
df_Path.Out.WP <- bind_rows(list(cbind(AllPathDR[[1]]$WP@result, "Timepoint" = rep("4h", dim(AllPathDR[[1]]$WP@result)[1])),
                                 cbind(AllPathDR[[2]]$WP@result, "Timepoint" = rep("8h", dim(AllPathDR[[2]]$WP@result)[1])),
                                 cbind(AllPathDR[[3]]$WP@result, "Timepoint" = rep("16h", dim(AllPathDR[[3]]$WP@result)[1])),
                                 cbind(AllPathDR[[4]]$WP@result, "Timepoint" = rep("24h", dim(AllPathDR[[4]]$WP@result)[1])),
                                 cbind(AllPathDR[[5]]$WP@result, "Timepoint" = rep("48h", dim(AllPathDR[[5]]$WP@result)[1])),
                                 cbind(AllPathDR[[6]]$WP@result, "Timepoint" = rep("72h", dim(AllPathDR[[6]]$WP@result)[1]))
))
df_Path.Out.WP$Timepoint <- as.factor(df_Path.Out.WP$Timepoint)


### BMD pathway output with the median BMD
entrez_id <- manif[match(df_All.Out.DRomic$id, manif$Gene_Symbol_TS),"Entrez.ID"]
gene_symbol <- manif[match(df_All.Out.DRomic$id, manif$Gene_Symbol_TS),"gene_symbol"]

df_All.Out.DRomic <- cbind(df_All.Out.DRomic, gene_symbol, entrez_id)
# write.table(df_All.Out.DRomic, "All_FilteredDRomic.txt", row.names = FALSE, sep = "\t")


## find MEDIAN by Pathway output and add value as column
findMedian <- function(x, BMDfile){
  
  Median <- list()
  
  for(i in 1:dim(x)[1]){
    
    lt <- strsplit(x[i, "geneID"], split = "/")
    tim <- x[i, 'Timepoint']
    df <- BMDfile[which(BMDfile$timepoint == tim), ]
    s <- df$BMD.zSD[match(lt[[1]], df$entrez_id)]
    Median[[i]] <- median(s)
    
  }
  
  unlist(Median)
}

# REACTOME
MedianReact <- findMedian(df_Path.Out.React, df_All.Out.DRomic)
df_Path.Out.React <- cbind(df_Path.Out.React, MedianReact)

# KEGG
MedianKEGG <- findMedian(df_Path.Out.KEGG, df_All.Out.DRomic)
df_Path.Out.KEGG <- cbind(df_Path.Out.KEGG, MedianKEGG)

# HALLMARK
MedianHALL <- findMedian(df_Path.Out.HALL, df_All.Out.DRomic)
df_Path.Out.HALL <- cbind(df_Path.Out.HALL, MedianHALL)

# WP
MedianWP <- findMedian(df_Path.Out.WP, df_All.Out.DRomic)
df_Path.Out.WP <- cbind(df_Path.Out.WP, MedianWP)


### Save Outputs
write.table(df_Path.Out.React, "All_ReactOut.txt", row.names = FALSE, sep = "\t")
write.table(df_Path.Out.KEGG, "All_KEGGOut.txt", row.names = FALSE, sep = "\t")
write.table(df_Path.Out.HALL, "All_HALLOut.txt", row.names = FALSE, sep = "\t")
write.table(df_Path.Out.WP, "All_WPOut.txt", row.names = FALSE, sep = "\t")






### For DEGs: ###
# while loop for res.* objects to create a list of data frames_________________________________________ START
i = 1
temp <- objects(pattern = "res\\.CS")
DEG.framesL <- list()

while(i <= length(temp)){
  
  cat(paste("Processing", i, "of", length(temp), "\n", sep = " "))
  s <- get(temp[i])
  df <- data.frame(s)
  gene_symbol <- manif[match(rownames(df), manif$Probe.Name),"gene_symbol"]
  entrez_id <- manif[match(rownames(df), manif$Probe.Name),"Entrez.ID"]
  
  df <- cbind(df, gene_symbol, entrez_id)
  DEG.framesL[[i]] <- df
  
  i = i + 1
}
names(DEG.framesL) <- temp

FindAGene <- sapply(DEG.framesL, function(df) grep("HAVCR1", df$gene_symbol))

# while loop for res.* objects to create a list of data frames_________________________________________ END
## Prepare MeSH hub, data frames. match manif(ests) with results tables, etc_________________________________________________________________________________________ END

# while loop for Gene Set Enrichment Analysis on the DEG.framesL_________________________________________ START
i = 1
temp <- length(DEG.framesL)
AllPath <- list()

while(i <= temp){
  
  cat(paste("Processing", i, "of", length(DEG.framesL), "\n", sep = " "))
  s <- DEG.framesL[[i]]
  
  sdeg <- s$entrez_id[which(s$padj <= 0.05)]    # for over-representation analysis (enrichXXX)
  
  geneLista <- s$log2FoldChange[which(s$padj <= 0.05)]
  names(geneLista) <- sdeg
  geneLista <- sort(geneLista, decreasing = TRUE)   # Gene Set Enrichment Analysis (gseXXX)
  
  REACT.gse <- gsePathway(geneLista)
  DO.gse <- gseDO(geneLista)
  KEGG.gse <- gseKEGG(geneLista)
  MeSH.gse <- gseMeSH(geneLista, MeSHDb = db, database='gendoo', category = 'C', maxGSSize = length(sdeg)-1)
  WP.gse <- gseWP(geneLista, organism = "Homo sapiens")
  
  AllPath[[i]] <- list(REA = REACT.gse,
                       DO = DO.gse,
                       KEGG = KEGG.gse,
                       MESH = MeSH.gse,
                       WP = WP.gse)
  
  i = i + 1
}
names(AllPath) <- gsub("res", "path", names(DEG.framesL))

# List of geneLists
i = 1
temp <- length(DEG.framesL)
AllLista <- list()

while(i <= temp){
  
  cat(paste("Processing", i, "of", length(DEG.framesL), "\n", sep = " "))
  s <- DEG.framesL[[i]]
  
  sdeg <- s$entrez_id[which(s$padj <= 0.05)]    # for over-representation analysis (enrichXXX)
  
  geneLista <- s$log2FoldChange[which(s$padj <= 0.05)]
  names(geneLista) <- sdeg
  geneLista <- sort(geneLista, decreasing = TRUE)   # Gene Set Enrichment Analysis (gseXXX)
  
  AllLista[[i]] <- geneLista
  
  i = i + 1
}
names(AllLista) <- names(DEG.framesL)

# List of plots
i = 1
temp <- length(DEG.framesL)
AllPlots <- list()

while(i <= temp){
  
  cat(paste("Processing", i, "of", length(DEG.framesL), "\n", sep = " "))
  s <- DEG.framesL[[i]]
  
  edox <- setReadable(AllPath[[i]]$REA, 'org.Hs.eg.db', 'ENTREZID')
  if(nrow(edox)>0){
    dplotR <- dotplot(edox, showCategory = 15, x = "GeneRatio", color = "p.adjust")
    hplotR <- heatplot(edox, foldChange = AllLista[[i]] , showCategory = 10)
    cplotR <- cnetplot(edox, foldChange = AllLista[[i]], showCategory = 20)
  } else {dplotR <- NULL
  hplotR <- NULL
  cplotR <- NULL}
  
  edox <- setReadable(AllPath[[i]]$DO, 'org.Hs.eg.db', 'ENTREZID')
  if(nrow(edox)>0){
    dplotD <- dotplot(edox, showCategory = 15, x = "GeneRatio", color = "p.adjust")
    hplotD <- heatplot(edox, foldChange = AllLista[[i]] , showCategory = 10)
    cplotD <- cnetplot(edox, foldChange = AllLista[[i]], showCategory = 20)
  } else {dplotD <- NULL
  hplotD <- NULL
  cplotD <- NULL}
  
  edox <- setReadable(AllPath[[i]]$KEGG, 'org.Hs.eg.db', 'ENTREZID')
  if(nrow(edox)>0){  
    dplotK <- dotplot(edox, showCategory = 15, x = "GeneRatio", color = "p.adjust")
    hplotK <- heatplot(edox, foldChange = AllLista[[i]] , showCategory = 10)
    cplotK <- cnetplot(edox, foldChange = AllLista[[i]], showCategory = 20)
  } else {dplotK <- NULL
  hplotK <- NULL
  cplotK <- NULL}
  
  edox <- setReadable(AllPath[[i]]$MESH, 'org.Hs.eg.db', 'ENTREZID')
  if(nrow(edox)>0){
    dplotM <- dotplot(edox, showCategory = 15, x = "GeneRatio", color = "p.adjust")
    hplotM <- heatplot(edox, foldChange = AllLista[[i]] , showCategory = 10)
    cplotM <- cnetplot(edox, foldChange = AllLista[[i]], showCategory = 20)
  } else {dplotM <- NULL
  hplotM <- NULL
  cplotM <- NULL}
  
  edox <- setReadable(AllPath[[i]]$WP, 'org.Hs.eg.db', 'ENTREZID')
  if(nrow(edox)>0){
    dplotW <- dotplot(edox, showCategory = 15, x = "GeneRatio", color = "p.adjust")
    hplotW <- heatplot(edox, foldChange = AllLista[[i]] , showCategory = 10)
    cplotW <- cnetplot(edox, foldChange = AllLista[[i]], showCategory = 20)
  } else {dplotW <- NULL
  hplotW <- NULL
  cplotW <- NULL}
  
  AllPlots[[i]] <- list(PlotsREA = list(dotREA = dplotR, heatREA = hplotR, netREA = cplotR),
                        PlotsDO = list(dotDO = dplotD, heatDO = hplotD, netDO = cplotD),
                        PlotsKEGG = list(dotKEGG = dplotK, heatKEGG = hplotK, netKEGG = cplotK),
                        PlotsMESH = list(dotMESH = dplotM, heatMESH = hplotM, netMESH = cplotM),
                        PlotsWP = list(dotWP = dplotW, heatWP = hplotW, netWP = cplotW))
  
  i = i + 1
}
names(AllPlots) <- names(DEG.framesL)


# while loop forGene Set Enrichment Analysis on the DEG.framesL____________________________________________ END


# List of significant entrez IDs in geneLists
i = 1
temp <- length(DEG.framesL)
AllNames <- list()

while(i <= temp){
  
  cat(paste("Processing", i, "of", length(DEG.framesL), "\n", sep = " "))
  s <- DEG.framesL[[i]]
  
  sdeg <- s$entrez_id[which(s$padj <= 0.05)]    # for over-representation analysis (enrichXXX)
  
  AllNames[[i]] <- sdeg
  
  i = i + 1
}
names(AllNames) <- names(DEG.framesL)


### Find a String in a list_______________________________________________START
library(tidyverse)
detect_string <- function(your_list, vector_strings){
  lapply(your_list, function(x) {
    if(TRUE %in% str_detect(x, paste(vector_strings, collapse = "|"))){
      TRUE
    } else {FALSE}
  })
}

