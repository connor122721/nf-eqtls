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

# Standardizing using MungeSumstats
if (args$liftover == TRUE) {
  
    print("Processing LiftOver!")
    param_lift="GRCh37"
    refor <- format_sumstats(path = gwas, 
                             ref_genome = param_lift, 
                             convert_ref_genome = "GRCh38", 
                             nThread = 8);
} else {
    param_lift="GRCh38"
    refor <- format_sumstats(path = gwas, 
                             ref_genome = param_lift,
                             nThread = 8)
}
                                  
print("Finished standardizing columns")

# Filter for significance and duplicated annotations
refor <- fread(refor)
refor.t <- data.table(refor[!duplicated(CHRBP_B38)])

# Restrict to all SNPs
gwas_SNPS <- unique(refor$CHRBP_B38)

# Split GWAS by chromosome and write 
foreach(i=1:length(uniqChrom)) %do% {
  print(uniqChrom[i])
  saveRDS(refor.t[CHROMOSOME==uniqChrom[i]], 
          file = paste("processed_", output_pre, "_", uniqChrom[i], ".rds", sep=""))
}

# Output entire liftover GWAS
# saveRDS(gwas, file = "processed_levin22_gwas_HF.rds")
