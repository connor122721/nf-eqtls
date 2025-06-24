#!/usr/bin/env Rscript
# Connor Murray
# Started 10.21.2024; edited 6.10.2025
# Reformat Splicing Data for TensorQTL

# Libraries
library(tidyverse)
library(data.table)
library(argparse)
library(foreach)

# Argument parser
parser <- ArgumentParser()
parser$add_argument("--metadata", required=TRUE, help="Path to metadata file")
parser$add_argument("--pca", required=TRUE, help="Path to splicing-PCA matrix")
parser$add_argument("--pca_snp", required=TRUE, help="Path to DNA-PCA matrix for ancestry")
args <- parser$parse_args()

#setwd("/standard/vol185/cphg_Manichaikul/users/csm6hg/");meta=fread("metadata/metadata_10_17_2024_CSM.txt");pca=fread("nextflow_dna/output/splicing/cluster/topchef.leafcutter.PCs.txt");dna=data.table(fread("nextflow_dna/output/pca/filt_dna_pc1_5_topchef.txt", header=T));norm.dt = fread("nextflow_dna/output/splicing/cluster/topchef.leafcutter.bed.gz")

# Metadata
meta <- data.table(fread(args$metadata, header = T))
print(args$metadata)
print(args$pca)
print(args$pca_snp)
# Read in PCA data
pca <- data.table(fread(args$pca, header = T))

# Make PCA friendly metadata
meta.pca <- data.table(t(meta[!SAMPLE_ID_TOR=="TOR238072"] %>% 
                        mutate(Gender=case_when(Gender == "Female" ~ 0,
                                                Gender == "Male" ~ 1,
                                                Gender == "" ~ 2),
                               Affected_NF=case_when(Affected_NF=="Affected" ~ 0,
                                                     Affected_NF=="Nonfailing" ~ 1),
                               CM_status=case_when(diagnosis_simple=="IDCM" ~ 0,
                                                   diagnosis_simple=="ICM" ~ 1,
                                                   TRUE ~ 2)) %>%  
                         select(sample=SAMPLE_ID_TOR, SAMPLE_ID_NWD, 
                                Gender, Age_at_collection, diagnosis_simple, Affected_NF)))
colnames(meta.pca) <- as.character(meta.pca[1,])

# Add ID column to match
meta.pca$ID = c("sample", "SAMPLE_ID_NWD", "Gender", 
                        "Age_at_collection", "diagnosis_simple", "Affected_NF")

# Include DNA-pca data
dna <- data.table(fread(args$pca_snp))
colnames(dna)[1] <- c("ID")

# Final Splicing-PCA; omit any columns with NAs
fin <- pca
t <- rbind(pca %>% mutate_if(is.double, as.character), 
           meta.pca, fill=T) %>% 
  select_if(~ !any(is.na(.)))

# Rename Splicing PCs
t$ID[t$ID %in% paste("PC", seq(1, 100), sep="")] <- paste("splice_PC", seq(1, 100), sep="")

# Replace RNA id with DNA identifier and merge w/ DNA PCA
colnames(t)[-1] <- as.character(t[ID=="SAMPLE_ID_NWD"])[-1]
t <- rbind(t %>% filter(!ID %in% c("sample", "SAMPLE_ID_NWD")), dna, fill=T) %>% 
  select_if(~ !any(is.na(.)))

###### Make format for tensorQTL ######

print(head(t))

# Output RNA PCs 1-100 for QTL saturation
foreach(i=1:100) %do% {
  
  # Message: i=10
  print(i)
  
  # List of PC dimensions
  PClist <- paste("splice_PC", seq(1, i), sep="")
  others <- t$ID[!t$ID %like% "splice"]
  rows = unique(c(PClist, others))
  
  # Remove samples with no age / sex # Affected group
  fin.dt <- na.omit(t[!ID=="diagnosis_simple"])
  colnames(fin.dt)[1] <- "x"
  
  # Restrict to specific combination of PCs for saturation assay
  fin1 <- fin.dt[x %in% rows]
  
  # Write covariates object
  write.table(x = fin1, 
              file = paste("topchef_cov_splicepc1_", i, "_06.23.25.txt", sep=""), 
              sep = "\t", 
              quote = F, 
              row.names = F)
  
}

#### Output RNAseq expression data


# Output final sample list
new.names = unique(colnames(t)[-1])
write.table(x = new.names, 
            file = "topchef_samples_06.24.25.txt", 
            sep = "\t", quote = F, row.names = F, col.names = F)
