# Connor Murray
# Started 12.10.2024; modifed 1.30.2025
# Make input file for LD pipeline
# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

# Libraries
library(data.table)
library(tidyverse)
library(foreach)

### Datasets & Setup ###
setwd("/standard/vol185/cphg_Manichaikul/users/csm6hg/")

gtf <- data.table(readRDS("genome_files/gencode.v34.GRCh38.ERCC.genes.collapsed.streamlined.RDS")) %>% 
  select(-c(V7:V9, file), gene_length=len)

# Read in cand gene list
cand <- fread("data/ClinGen_DCM_ARVC_GeneIDs.txt") %>% 
  left_join(gtf, by=c("common_gene")) %>% 
  select(c(chrom, common_gene, start, stop, gene_edit, class, evidence))

cand2 <- fread("nextflow_dna/output/coloc/coloc_eqtl_candidates_full.txt") %>% 
  filter(PP.H4 >= 0.5) %>% 
  left_join(gtf, by=c("common_gene", "start", "stop", "chrom")) %>% 
  select(chrom, common_gene, start, stop, gene_edit)

# Rbind
cand.fin <- data.table(rbind(cand, cand2, fill=T))

# Make LD parameter list
dist = 1e6
cand.dt <- cand.fin %>% 
  mutate(startLD = start-dist,
         stopLD = stop+dist)

# Output
fin <- cand.dt %>% select(gene_edit, common_gene, chrom, startLD, stopLD, evidence)
write_delim(fin, file = "data/candgenes.ldlist.txt", delim = "\t")
