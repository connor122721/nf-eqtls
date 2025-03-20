# Connor Murray
# Started 12.3.2024; edited 3.20.2025
# Standardize GWAS dataset and LiftOver from Hg19 -> Hg38
# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

# Libraries
library(data.table)
library(tidyverse)
library(foreach)
library(MungeSumstats)
library(argparse)

# Argument parser
parser <- ArgumentParser()
parser$add_argument("--gwas", required=TRUE, help="GWAS summary statistics.")
parser$add_argument("--prefix", required=TRUE, help="Prefix of output files.")
parser$add_argument("--liftover", help="Statment to liftover from hg19 -> hg38.")
args <- parser$parse_args()

gwas_input <- fread(args$gwas, header=T)
output_pre <- args$prefix
# gwas_input=fread("../data/HF-multiancestry-maf0.01.tsv.gz", header=T)
# gwas_input=fread("/standard/vol185/cphg_Manichaikul/users/csm6hg/data/DCM_GWAS/Jurgens_DCM_GWAS_META_BiobanksOnly.tsv.gz"); liftover="FALSE"

### Datasets & Setup ###

# Make standardized column names
newcols = data.table(cols=colnames(gwas_input)) %>% 
  mutate(new=case_when(cols=="base_pair_location"~"BP",
                       cols=="pos"~"BP",
                       cols=="POS"~"BP",
                       cols=="effect_allele_frequency"~"af",
                       cols=="freq"~"af",
                       cols=="EAFREQ"~"af",
                       cols=="chromosome"~"CHR",
                       cols=="p"~"p_value",
                       cols=="P"~"p_value",
                       cols=="SE"~"standard_error",
                       cols=="se"~"standard_error",
                       cols=="b"~"beta",
                       cols=="BETA"~"beta",
                       cols=="effect_allele"~"A1",
                       cols=="EA"~"A1",
                       cols=="NEA"~"A2",
                       cols=="other_allele"~"A2",
                       cols=="rsID"~"SNP",
                       cols=="variant_id"~"SNP",
                       TRUE ~ cols))

colnames(gwas_input) <- newcols$new

# Read in GWAS 
gwas <- data.table(gwas_input %>% 
                     mutate(snpID=paste("chr", CHR, ":", BP, sep=""),
                            maf=case_when(af > 0.5 ~ 1 - af,
                                          TRUE ~ af),
                            chromosome=paste("chr", CHR, sep="")))

print("Finished standardizing columns")

# Check if liftover should be processed
if (args$liftover == TRUE) {
  
  # LiftOver
  library(rtracklayer)
  print("Processing LiftOver!")
  
  # Download liftover files
  # wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
  # wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
  
  # Load liftOver chain file
  chain <- import.chain("/standard/vol185/cphg_Manichaikul/users/csm6hg/genome_files/hg19ToHg38.over.chain")
  
  # Convert to GRanges
  gr <- GRanges(seqnames = gwas$CHR,
                ranges = IRanges(start = gwas$BP, 
                                 end = gwas$BP))
  
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
  uniqChrom = na.omit(unique(gwas$chr_hg38))
} else {
  
  print("Skipping LiftOver!")
  
  gwas[, chr_hg38 := CHR]
  gwas[, pos_hg38 := BP]
  gwas[, snpID_hg38 := snpID]
  uniqChrom = na.omit(unique(gwas$chromosome))
}

# Check Effect Allele vs. Reference Allele
library(BSgenome.Hsapiens.UCSC.hg38)
ref_genome <- BSgenome.Hsapiens.UCSC.hg38

#gwasTest=gwas[1:1000]

# For each SNP, extract the reference base from the genome using chr_hg38 and pos_hg38
gwasTest[, ref_allele := mapply(function(chr, pos) {
  if (chr %in% seqnames(ref_genome)) {
    as.character(subseq(ref_genome[[chr]], start = pos, end = pos))
  } else {
    NA_character_
  }
}, chromosome, pos_hg38)]

# Compare the reported effect allele (A1) to the reference allele
gwas[, allele_flip := ifelse(toupper(A1) != toupper(ref_allele) & toupper(A2) == toupper(ref_allele), TRUE, FALSE)]

if ("beta" %in% colnames(gwasTest)) {
  gwas[allele_flip == TRUE, c("A1", "A2") := list(A2, A1)]
  gwas[allele_flip == TRUE, beta := -beta]
} else {
  gwas[allele_flip == TRUE, c("A1", "A2") := list(A2, A1)]
}

# How many SNPs match the reference?
print("Summary of effect allele check (TRUE means A1 matches the reference allele):")
print(table(gwas$effect_allele_correct))

# Filter for significance and duplicated annotations
gwas.t <- data.table(gwas[!duplicated(snpID_hg38)])

# Restrict to all SNPs
gwas_SNPS <- unique(gwas$snpID_hg38)

# Split GWAS by chromosome and write 
foreach(i=1:length(uniqChrom)) %do% {
  print(uniqChrom[i])
  saveRDS(gwas.t[chromosome == uniqChrom[i]], 
          file = paste("processed_", output_pre, "_", uniqChrom[i], ".rds", sep=""))
}

# Output entire liftover GWAS
# saveRDS(gwas, file = "processed_levin22_gwas_HF.rds")