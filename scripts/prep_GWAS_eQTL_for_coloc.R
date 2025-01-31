# Connor Murray
# Started 12.3.2024; edited 1.27.2025
# Standardize GWAS dataset and LiftOver from Hg19 -> Hg38
# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

# Libraries
library(data.table)
library(tidyverse)
library(foreach)
library(argparse)

# Argument parser
parser <- ArgumentParser()
parser$add_argument("--gwas", required=TRUE, help="GWAS summary statistics.")
parser$add_argument("--prefix", required=TRUE, help="Prefix of output files.")
parser$add_argument("--skip_liftover", help="Statment to skip liftover from hg19 -> hg38.")
args <- parser$parse_args()

gwas_input <- fread(args$gwas, header=T)
output_pre <- args$prefix
# gwas_input=fread("../data/HF-multiancestry-maf0.01.tsv.gz", header=T)
# gwas_input=fread("../../../../data/DCM_GWAS/Jurgens_DCM_GWAS_META_BiobanksOnly.tsv.gz")

### Datasets & Setup ###

# Make standardized column names
newcols = data.table(cols=colnames(gwas_input)) %>% 
  mutate(new=case_when(cols=="base_pair_location"~"pos",
                       cols=="BP"~"pos",
                       cols=="POS"~"pos",
                       cols=="effect_allele_frequency"~"af",
                       cols=="freq"~"af",
                       cols=="EAFREQ"~"af",
                       cols=="CHR"~"chromosome",
                       cols=="p"~"p_value",
                       cols=="P"~"p_value",
                       cols=="SE"~"standard_error",
                       cols=="se"~"standard_error",
                       cols=="b"~"beta",
                       cols=="BETA"~"beta",
                       TRUE ~ cols))

colnames(gwas_input) <- newcols$new

# Read in GWAS 
gwas <- data.table(gwas_input %>% 
              mutate(snpID=paste("chr", chromosome, ":", pos, sep=""),
                     maf=case_when(af > 0.5 ~ 1 - af,
                                   TRUE ~ af)))
                                  
print("Finished standardizing columns")

# Check if liftover should be skipped
if (is.null(args$skip_liftover)) {

  # LiftOver
  library(rtracklayer)

  # Download liftover files
  # wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
  # wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

  # Load liftOver chain file
  chain <- import.chain("/standard/vol185/cphg_Manichaikul/users/csm6hg/genome_files/hg19ToHg38.over.chain")

  # Convert to GRanges
  gr <- GRanges(seqnames = paste("chr", gwas$chromosome, sep=""),
                ranges = IRanges(start = gwas$pos, 
                                 end = gwas$pos))

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
} else {
  gwas[, chr_hg38 := chromosome]
  gwas[, pos_hg38 := pos]
  gwas[, snpID_hg38 := snpID]
  uniqChrom=na.omit(unique(gwas$chromosome))
}

# Filter for significance and duplicated annotations
gwas.t <- data.table(gwas[!duplicated(snpID_hg38)])

# Restrict to all SNPs
gwas_SNPS <- unique(gwas$snpID_hg38)

# Split GWAS by chromosome and write 
foreach(i=1:length(uniqChrom)) %do% {
  print(uniqChrom[i])
  saveRDS(gwas.t[chr_hg38==uniqChrom[i]], 
          file = paste("processed_", output_pre, "_", uniqChrom[i], ".rds", sep=""))
}

# Output entire liftover GWAS
# saveRDS(gwas, file = "processed_levin22_gwas_HF.rds")
