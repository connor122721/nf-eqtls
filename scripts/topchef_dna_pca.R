#!/usr/bin/env Rscript
# Performs PCA on GDS file and generates PCA plots and TensorQTL input

# Libraries
library(data.table)
library(tidyverse)
library(SeqArray)
library(SNPRelate)
library(ggforce)
library(patchwork)
library(parallel)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check for correct number of arguments
if(length(args) != 7){
  stop("Usage: perform_pca.R <input_gds> <metadata_file> <related_individuals> <pca_rds> <pca_plot> <tensorqtl_pca> <outliers_nam>")
}

input_gds <- args[1]
metadata_file <- args[2]
related_individuals <- args[3]
pca_rds <- args[4]
pca_plot <- args[5]
tensorqtl_pca <- args[6]
outlier_output <- args[7]

# Register cores
threads <- parallel::detectCores()
seqParallelSetup(threads)

# Open GDS
genofile <- seqOpen(input_gds)
seqResetFilter(genofile)

# Samples full
samps <- seqGetData(genofile, var.name = "sample.id")

# Sample metadata
fin <- data.table(fread(metadata_file))
related <- data.table(fread(related_individuals, header = FALSE))
fin <- fin[!SAMPLE_ID_NWD %in% related$V1]
seqSetFilter(genofile, sample.id = fin$SAMPLE_ID_NWD)

# SNPs metadata
dt <- data.table(
  variant.id = seqGetData(genofile, var.name = "variant.id"),
  position = seqGetData(genofile, var.name = "position"),
  chrom = seqGetData(genofile, var.name = "chromosome"),
  ac = seqAlleleCount(genofile, ref.allele=0L)) %>% 
  mutate(
    af = ac / (length(samps) * 2),
    maf = case_when(
      af > 0.5 ~ 1 - af,
      TRUE ~ af))

# Run PCA
ccm_pca <- snpgdsPCA(
  genofile, 
  autosome.only = FALSE, 
  sample.id = as.character(fin$SAMPLE_ID_NWD),
  maf = 0.05, 
  missing.rate = 0,
  remove.monosnp = TRUE,
  num.thread = threads)

# Output PCA
saveRDS(ccm_pca, file = pca_rds)

# Scree plot data
scree_data <- data.frame(
  var = c(ccm_pca$varprop)[1:15] * 100,
  id = 1:15)

# Merge by Sample
pca <- data.table(
  sample = ccm_pca$sample.id, 
  PC1 = ccm_pca$eigenvect[,1],
  PC2 = ccm_pca$eigenvect[,2],
  PC3 = ccm_pca$eigenvect[,3],
  PC4 = ccm_pca$eigenvect[,4],
  PC5 = ccm_pca$eigenvect[,5],
  PC6 = ccm_pca$eigenvect[,6],
  PC7 = ccm_pca$eigenvect[,7],
  PC8 = ccm_pca$eigenvect[,8],
  PC9 = ccm_pca$eigenvect[,9],
  PC10 = ccm_pca$eigenvect[,10])

# Remove missing samples
pca <- data.table(pca[sample %in% fin$SAMPLE_ID_NWD])

# Merge PCA and metadata
pca <- data.table(merge(pca, fin, by.x = "sample", by.y = "SAMPLE_ID_NWD"))

######## Outlier detection and filtering ########

# Mahalanobis distance calculation
mal <- mahalanobis(
  x = pca %>% select(PC1:PC5),
  center = colMeans(pca %>% select(PC1:PC5)),
  cov = cov(pca %>% select(PC1:PC5)))

# Create new column in data frame to hold p-value for each Mahalanobis distance
p <- pchisq(mal, df = 4, lower.tail = FALSE)

# Mal distance data - adjust p-values
mal.dt <- data.table(
  pca %>% select(sample, Affected_NF), 
  mahal = mal, 
  pval = p,
  padj = p.adjust(p, method = "BH"))

outliers <- mal.dt[padj < 0.05]$sample

# Common theme element
themei <- {
  theme_bw() + 
    theme(
      axis.title.x = element_text(face = "bold", size = 16),
      axis.text.x = element_text(face = "bold", size = 14),
      axis.title.y = element_text(face = "bold", size = 16),
      axis.text.y = element_text(face = "bold", size = 14),
      legend.text = element_text(face = "bold", size = 14),
      legend.title = element_text(face = "bold", size = 16))}

# Fix self-reported ancestry
pca <- data.table(pca %>% 
                   mutate(
                     RaceSimp = case_when(
                       Race %in% c("", "?", "N/A", "Unknown or not reported", "Other") ~ "Unknown",
                       Race %in% c("American Indian/Alaska Native", "Native Hawaiian or other Pacific Islander") ~ "Native Am.",
                       Race == "Black or African American" ~ "African Am.",
                       Race == "More than one race" ~ "Admixed",
                       TRUE ~ Race)))

# PC 1/2 Plot
pc12 <- {
  pca[!sample %in% outliers] %>%
    ggplot(aes(x = PC1, y = PC2, color = RaceSimp)) +
    geom_point(size = 6, alpha = 0.6) +
    theme_minimal() + 
    labs(
      x = paste("PC1 (", round(ccm_pca$varprop[[1]], digits = 3) * 100, " %)", sep = ""),
      y = paste("PC2 (", round(ccm_pca$varprop[[2]], digits = 3) * 100, " %)", sep = "")) +
    theme(
      strip.text = element_text(face = "bold.italic", size = 16),
      legend.text = element_text(size = 16, face = "bold.italic"),
      legend.title = element_text(face = "bold", size = 18),
      axis.text.x = element_text(face = "bold", size = 18),
      axis.text.y = element_text(face = "bold", size = 18),
      axis.title.x = element_text(face = "bold", size = 20),
      axis.title.y = element_text(face = "bold", size = 20),
      axis.title = element_text(face = "bold", size = 20))
}

# PC 2/3 Plot
pc23 <- {
  pca[!sample %in% outliers] %>%
    ggplot(aes(x = PC2, y = PC3, color = RaceSimp)) +
    geom_point(size = 6, alpha = 0.6) +
    theme_minimal() + 
    labs(
      x = paste("PC2 (", round(ccm_pca$varprop[[2]], digits = 3) * 100, " %)", sep = ""),
      y = paste("PC3 (", round(ccm_pca$varprop[[3]], digits = 3) * 100, " %)", sep = "")) +
    theme(
      strip.text = element_text(face = "bold.italic", size = 16),
      legend.text = element_text(size = 16, face = "bold.italic"),
      legend.title = element_text(face = "bold", size = 18),
      axis.text.x = element_text(face = "bold", size = 18),
      axis.text.y = element_text(face = "bold", size = 18),
      axis.title.x = element_text(face = "bold", size = 20),
      axis.title.y = element_text(face = "bold", size = 20),
      axis.title = element_text(face = "bold", size = 20))
}

# PC 3/4 Plot
pc34 <- {
  pca[!sample %in% outliers] %>%
    ggplot(aes(x = PC3, y = PC4, color = RaceSimp)) +
    geom_point(size = 6, alpha = 0.6) +
    theme_minimal() + 
    labs(
      x = paste("PC3 (", round(ccm_pca$varprop[[3]], digits = 3) * 100, " %)", sep = ""),
      y = paste("PC4 (", round(ccm_pca$varprop[[4]], digits = 3) * 100, " %)", sep = "")) +
    theme(
      strip.text = element_text(face = "bold.italic", size = 16),
      legend.text = element_text(size = 16, face = "bold.italic"),
      legend.title = element_text(face = "bold", size = 18),
      axis.text.x = element_text(face = "bold", size = 18),
      axis.text.y = element_text(face = "bold", size = 18),
      axis.title.x = element_text(face = "bold", size = 20),
      axis.title.y = element_text(face = "bold", size = 20),
      axis.title = element_text(face = "bold", size = 20))
}

# PC 4/5 Plot
pc45 <- {
  pca[!sample %in% outliers] %>%
    ggplot(aes(x = PC4, y = PC5, color = RaceSimp)) +
    geom_point(size = 6, alpha = 0.6) +
    theme_minimal() + 
    labs(
      x = paste("PC4 (", round(ccm_pca$varprop[[4]], digits = 3) * 100, " %)", sep = ""),
      y = paste("PC5 (", round(ccm_pca$varprop[[5]], digits = 3) * 100, " %)", sep = "")) +
    theme(
      strip.text = element_text(face = "bold.italic", size = 16),
      legend.text = element_text(size = 16, face = "bold.italic"),
      legend.title = element_text(face = "bold", size = 18),
      axis.text.x = element_text(face = "bold", size = 18),
      axis.text.y = element_text(face = "bold", size = 18),
      axis.title.x = element_text(face = "bold", size = 20),
      axis.title.y = element_text(face = "bold", size = 20),
      axis.title = element_text(face = "bold", size = 20))
}

# PC 5/6 Plot
pc56 <- {
  pca[!sample %in% outliers] %>%
    ggplot(aes(x = PC5, y = PC6, color = RaceSimp)) +
    geom_point(size = 6, alpha = 0.6) +
    theme_minimal() + 
    labs(
      x = paste("PC5 (", round(ccm_pca$varprop[[5]], digits = 3) * 100, " %)", sep = ""),
      y = paste("PC6 (", round(ccm_pca$varprop[[6]], digits = 3) * 100, " %)", sep = "")) +
    theme(
      strip.text = element_text(face = "bold.italic", size = 16),
      legend.text = element_text(size = 16, face = "bold.italic"),
      legend.title = element_text(face = "bold", size = 18),
      axis.text.x = element_text(face = "bold", size = 18),
      axis.text.y = element_text(face = "bold", size = 18),
      axis.title.x = element_text(face = "bold", size = 20),
      axis.title.y = element_text(face = "bold", size = 20),
      axis.title = element_text(face = "bold", size = 20))
}

# Scree plot
scree <- {
  ggplot(scree_data, aes(x = id, y = var)) +
    geom_point(size = 2) +
    geom_vline(xintercept = 5, linetype = 2) +
    geom_line() +
    labs(x = "PC Loading", y = "PC Variance %") +
    themei
}

# Composite figure
pc_fin <- (((pc12 | pc23 | pc34) / (pc45 | pc56 | scree)) + 
            plot_layout(guides = 'collect'))

# Save composite PCA plot
ggsave(plot = pc_fin, filename = pca_plot, width = 16, height = 8)

# Output PCA for TensorQTL
write.table(
  t(pca %>% select(sample, PC1:PC5)),
  file = tensorqtl_pca,
  quote = FALSE,
  row.names = TRUE,
  col.names = FALSE,
  sep = '\t')

# Output outlier individuals
write.table(
  outliers,
  file = outlier_output,
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  sep = '\t')
