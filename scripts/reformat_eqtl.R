#!/usr/bin/env Rscript
# Connor Murray
# Started 10.21.2024

# Libraries
library(tidyverse)
library(data.table)
library(foreach)
library(FactoMineR)
library(PCAForQTL)
library(argparse)

# Argument parser
parser <- ArgumentParser()
parser$add_argument("--metadata", required=TRUE, help="Path to metadata file")
parser$add_argument("--gene_gtf", required=TRUE, help="Path to GTF file")
parser$add_argument("--related_individuals", required=TRUE, help="Path to related individuals file")
parser$add_argument("--rna_outliers", required=TRUE, help="Path to RNA outliers file")
parser$add_argument("--dna_outliers", required=TRUE, help="Path to DNA outliers file")
parser$add_argument("--norm_tmm", required=TRUE, help="Path to normalized TMM RNAseq data")
parser$add_argument("--pca_tmm", required=TRUE, help="Path to PCA TMM data")
parser$add_argument("--pca_snp", required=TRUE, help="Path to SNP PCA data")
args <- parser$parse_args()

#setwd("/standard/vol185/cphg_Manichaikul/users/csm6hg/")
#metadata="metadata/metadata_10_17_2024_CSM.txt"
#gene_gtf="genome_files/gencode.v34.GRCh38.ERCC.genes.collapsed.TSS.bed"
#related_individuals="nextflow_dna/output_test_new/king/relatedIndividuals.txt"
#norm_tmm="nextflow_dna/output_test_new/rna/norm_tmm.tsv"
#pca_tmm="nextflow_dna/output_test_new/rna/pca_tmm.tsv"
#pca_snp="nextflow_dna/output_test_new/pca/filt_dna_pc1_5_topchef.txt"

# Metadata
meta <- data.table(fread(args$metadata, header = T))

# Read in normalized TMM RNAseq data
norm_tmm <- data.table(fread(args$norm_tmm, header = T))

# Read in PCA TMM data
pca_tmm <- data.table(fread(args$pca_tmm, header = T))

# Make PCA friendly metadata
meta.pca <- data.table(meta[!SAMPLE_ID_TOR=="TOR238072"] %>% 
                         select(sample=SAMPLE_ID_TOR, SAMPLE_ID_NWD, 
                                Gender, Age_at_collection, diagnosis_simple, Affected_NF))
rownames(meta.pca) <- meta.pca$sample

# Final RNA PCA
fin <- pca_tmm
colnames(fin)[1] <- "SAMPLE_ID_TOR"
pcrnacols = colnames(fin)[colnames(fin)%like%"PC"]
colnames(fin)[colnames(fin)%like%"PC"] = paste("RNA_", pcrnacols, sep="") 

# Add equivalent ID
fin <- na.omit(fin %>% 
  left_join(meta %>% 
              select("SAMPLE_ID_TOR", "SAMPLE_ID_NWD", "Gender")))

# Read in GTF
gtf_tss <- data.table(fread(args$gene_gtf))
gtf_tss[, gene_edit := sub("\\..*", "", gene_id)]
norm_tmm[, gene_edit := sub("\\..*", "", Name)]

# Close relatives in dataset (as inferred by king analysis) and outliers 
out.rels <- unlist(fread(args$related_individuals, header = F))
out.dna <- unlist(fread(args$dna_outliers, header = F))
out.rna <- unlist(fread(args$rna_outliers, header = F))

# Combine gene info with normalized RNA
filt_norm_counts_gene <- data.table(norm_tmm %>% 
                    left_join(gtf_tss %>% 
                                select(gene_edit, chr, start, end), 
                                      by=c("gene_edit")))
  

###### Make format for tensorQTL bed ######

#### DNA data
dna <- data.table(t(fread(args$pca_snp, header=F)))[-1]
colnames(dna) <- c("sample","PC1","PC2","PC3","PC4","PC5")

# Samples in Bed File

# All outliers
outs <- unique(c(out.rels, out.dna, out.rna))

# Hard code samples and remove outliers/related individuals
dna <- data.table(dna[!sample %in% outs])

# Output RNA PCs 1-100 for eQTL saturation
foreach(i=1:100) %do% {
  
  # Message: i=10
  print(i)
  
  # List of PC dimensions
  PClist <- paste("RNA_PC", seq(1, i), sep="")
  
  # RNA
  rnai <- data.frame(fin[SAMPLE_ID_NWD %in% dna$sample] %>% 
                  select(SAMPLE_ID_NWD, PClist, Gender, Affected_NF))
  
  # Merging - first round
  fin.dt <- data.table(rnai %>% 
                left_join(meta.pca %>% 
                            select(-c(sample)), 
                          by=c("SAMPLE_ID_NWD", "Gender", "Affected_NF")) %>% 
                left_join(dna, 
                          by=c("SAMPLE_ID_NWD" = "sample")) %>% 
                 mutate(Gender=case_when(Gender == "Female" ~ 0,
                                         Gender == "Male" ~ 1,
                                         Gender == "" ~ 2),
                        Affected_NF=case_when(Affected_NF=="Affected" ~ 0,
                                              Affected_NF=="Nonfailing" ~ 1),
                        CM_status=case_when(diagnosis_simple=="IDCM" ~ 0,
                                            diagnosis_simple=="ICM" ~ 1,
                                            TRUE ~ 2)) %>% 
                  select(-c(diagnosis_simple), sample=SAMPLE_ID_NWD))
  
  # Remove samples with no age / sex # Affected group
  fini <- t(na.omit(fin.dt))
  colnames(fini) <- fini[1,]
  fin1 <- data.table(fini[-1,])
  fin1$x <- rownames(fini[-1,])
  
  # Write covariates object
  write.table(x = fin1 %>% 
                select(c(x, contains("NWD"))), 
              file = paste("topchef_cov_RNApc1_", i, "_1.15.25.txt", sep=""), 
              sep = "\t", 
              quote = F, 
              row.names = F)
  
  dim(fin1)
  
}

# Switch sample names to those in the genetic data, also remove Y chromosome genes
out <- data.table(filt_norm_counts_gene[!chr %in% "chrY"] %>% 
                    select(chr, start, end, gene_edit, contains("TOR")))

out_s <- data.table(colnames(out)[-(1:4)]) %>% 
  left_join(meta %>% select(contains("SAMPLE_ID")), 
            by=c("V1"="SAMPLE_ID_TOR"))

colnames(out)[-(1:4)] <- out_s$SAMPLE_ID_NWD

sampi <- rownames(t(fin1))[rownames(t(fin1))%like%"NWD"]

#### Output RNAseq expression data
test <- out %>% 
  select("#chr"=chr, start, end, phenotype_id=gene_edit, contains(sampi))

write.table(test, 
            file = "filt_rnaseq_norm_topchef_unrelated.bed", 
            quote = F, row.names = F, sep="\t")

# Output final sample list
write.table(x = sampi, 
            file = "topchef_samples_1_15_25.txt", 
            sep = "\t", quote = F, row.names = F, col.names = F)

### Run Elbow Test ###

# Wide to long
pcmat <- t(filt_norm_counts_gene %>% select(-c(chr, start, end)))
colnames(pcmat) <- filt_norm_counts_gene$gene_edit

# Convert all elements to numeric
pcmat <- apply(pcmat, 2, as.numeric)

# Elbow test 
resultRunElbow <- PCAForQTL::runElbow(X = na.omit(pcmat))
print(resultRunElbow)

# Output elbow test results
write.table(x = resultRunElbow, 
            file = "rna_pca_elbow_best_k", 
            sep = "\t", quote = F, row.names = F, col.names = F)
