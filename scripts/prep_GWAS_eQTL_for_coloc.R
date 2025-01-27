# Connor Murray
# Started 12.3.2024; edited 1.27.2025
# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

# Libraries
library(data.table)
library(tidyverse)
library(foreach)
library(argparse)

# Argument parser
parser <- ArgumentParser()
parser$add_argument("--gwas", required=TRUE, help="GWAS summary statistics.")
args <- parser$parse_args()

gwas_input <- fread(args$gwas, header=T)
# gwas_input=fread("data/HF-multiancestry-maf0.01.tsv.gz", header=T)

### Datasets & Setup ###

# Read in GWAS 
gwas <- data.table(gwas_input %>% 
              mutate(snpID=paste("chr", chromosome, ":", base_pair_location, sep=""),
                     maf=case_when(effect_allele_frequency > 0.5 ~ 1 - effect_allele_frequency,
                                   TRUE ~ effect_allele_frequency)))

# LiftOver
library(rtracklayer)
# wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
# wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

# Load liftOver chain file
chain <- import.chain("/standard/vol185/cphg_Manichaikul/users/csm6hg/genome_files/hg19ToHg38.over.chain")

# Convert to GRanges
gr <- GRanges(seqnames = paste("chr", gwas$chromosome, sep=""),
              ranges = IRanges(start = gwas$base_pair_location, 
                               end = gwas$base_pair_location))

# Apply liftOver
lifted_gr <- liftOver(gr, chain)

# Handle unmapped positions by setting missing values to NA
lifted_gr_unlisted <- unlist(lifted_gr)  # Only mapped coordinates
mapped_indices <- which(lengths(lifted_gr) > 0)  # Indices of successfully mapped positions

# Create new columns with default NA values
gwas[, chr_hg38 := NA_character_]
gwas[, pos_hg38 := NA_integer_]

# Fill in mapped positions
gwas$chr_hg38[mapped_indices] <- as.character(seqnames(lifted_gr_unlisted))
gwas$pos_hg38[mapped_indices] <- start(lifted_gr_unlisted)
gwas[, snpID_hg38 := paste0(chr_hg38, ":", pos_hg38)]
uniqChrom=na.omit(unique(gwas$chr_hg38))

# Filter for significance and duplicated annotations
gwas.t <- data.table(gwas[!duplicated(snpID_hg38)])

# Restrict to all SNPs
gwas_SNPS <- unique(gwas$snpID_hg38)

# Split GWAS by chromosome and write 
foreach(i=1:length(uniqChrom)) %do% {
  print(uniqChrom[i])
  saveRDS(gwas.t[chr_hg38==uniqChrom[i]], 
          file = paste("processed_levin22_gwas_HF_", uniqChrom[i], ".rds", sep=""))
}

# Output
saveRDS(gwas, file = "processed_levin22_gwas_HF.rds")
