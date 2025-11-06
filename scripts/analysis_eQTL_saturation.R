# Connor Murray
# Started 10.21.2024; edit 9.11.2025
# Analyzing TOPchef saturation of significant e/sQTLs 

# Libraries
library(tidyverse)
library(stringi)
library(data.table)
library(patchwork)
library(future.apply)
library(argparse)

# Set up parallel backend, using all available cores
plan(multisession, workers = 4)

# Argument parser
parser <- ArgumentParser()
parser$add_argument("--list_of_eqtls", required=TRUE, help="List of TensorQTL output to analyze.")
parser$add_argument("--gtf", required=TRUE, help="Streamlined GTF object.")
parser$add_argument("--sqtl", action='store_true', required=FALSE, help="Make script usable with sQTL object.")
args <- parser$parse_args()

#setwd("/scratch/csm6hg/nextflow_dna/work/20/5029b6c2279e682ae5a7f7f68866f2"); qtl.files=fread("cis_sqtl.list", header=F)
#gene.gtf=read_rds("/standard/vol185/cphg_Manichaikul/users/csm6hg/genome_files/gencode.v34.GRCh38.ERCC.genes.collapsed.streamlined.RDS"); sqtl=TRUE

qtl.files <- fread(args$list_of_eqtls, header=F)
gtf <- args$gtf

# List reading files
qtl <- rbindlist(future_lapply(qtl.files$V1, function(file) {
  #data <- arrow::read_parquet(file)
  data <- fread(file)
  data[, .id := basename(file)]  # Add the file name as the .id column
  print(basename(file))
  return(data)
}))

#arrow::read_parquet(qtl.files[1])

# Add which run 
qtl.dt <- data.table(qtl %>% 
           mutate(maxPC = as.numeric(str_remove_all(str_remove_all(tstrsplit(.id, "_")[[3]], ".cis"), "MaxPC")),
                  chrom = as.factor(tstrsplit(.id, "_")[[2]])) %>% 
           select(-c(.id)))

# FDR pvalue adjusting
qtl.dt <- data.table(qtl.dt %>% 
           group_by(maxPC) %>% 
           mutate(padj = -log10(p.adjust(pval_nominal, method = "fdr"))))

# Read in GTF
gene.gtf <- data.table(readRDS(gtf))

# Conditional statement for gene metadata
if (args$sqtl == "TRUE") {
  
  # sQTLs add gene metadata
  print("sQTL Processing!")
  qtl.dt <- data.table(qtl.dt %>% 
             mutate(phenotype_id = tstrsplit(tstrsplit(phenotype_id, ":")[[5]], ".", fixed=T)[[1]]) %>%  
             left_join(gene.gtf %>% 
                      select(-c(file, V9)), 
                      by = c("phenotype_id" = "gene_edit", 
                             "chrom"="chrom")))
  prefix="sqtl"
  
  } else {
    
  # eQTLs add gene metadata
  print("eQTL Processing!")
  qtl.dt <- data.table(qtl.dt %>% 
             left_join(gene.gtf %>% 
                      select(-c(file, V9)), 
                      by = c("phenotype_id" = "gene_edit", 
                             "chrom"="chrom")))
  prefix="eqtl"
}

# Output
saveRDS(qtl.dt, file = paste0(prefix, ".rna.saturation.rds"))
# qtl.dt <- read_rds("sqtl.rna.saturation.rds")

#### Analysis ####

# Common theme element
themei <- {
  theme_bw() + 
    theme(axis.title.x = element_text(face = "bold", size = 16),
          axis.text.x = element_text(face = "bold", size = 14),
          axis.title.y = element_text(face = "bold", size = 16),
          axis.text.y = element_text(face = "bold", size = 14),
          legend.text = element_text(face = "bold", size = 14),
          legend.title = element_text(face = "bold", size = 16),
          aspect.ratio = 1) 
}

# Pvalue cutoff
cutoff = 0.05

### Saturation test

# Number of sig eQTLs and rate of change
qtl.rate <- data.table(qtl.dt %>% 
                 group_by(maxPC) %>% 
                 summarize(n = sum(pval_perm < cutoff)) %>% 
                 arrange(maxPC) %>%  
                 mutate(roc = (n - lag(n)) / (maxPC - lag(maxPC)),
                        roc_per = (n/lag(n))-1) %>%
                 mutate(roc2 = (roc - lag(roc)) / (maxPC - lag(maxPC))))

if (args$sqtl == "TRUE") {
  
  # find the index that maximizes second derivative in absolute value:
  best_k <- which.max(qtl.rate$n)
  best_maxPC <- qtl.rate$maxPC[best_k]
  best_maxPC
  
} else {

  # Correlation of significant genes
  qtl.cor <- data.table(qtl.dt %>%
    mutate(sig = case_when(pval_perm < cutoff ~ 1,
                           TRUE ~ 0)) %>%  # Filtering for significant genes
    select(maxPC, phenotype_id, sig) %>%  # Select relevant columns
    spread(key = maxPC, value = sig)) # Make long to wide matrix
  
  # Pairwise correlations
  cor_matrix <- cor(qtl.cor[,-1], use = "pairwise.complete.obs")
  upper_triangle <- cor_matrix
  upper_triangle[lower.tri(cor_matrix, diag = TRUE)] <- NA  # Mask lower triangle and diagonal
  
  # Convert to a data table for viewing
  cor_list <- as.data.table(as.table(cor_matrix))
  
  # find the index that maximizes second derivative in absolute value:
  best_k <- which.max(qtl.rate$n)
  best_maxPC <- qtl.rate$maxPC[best_k]
  best_maxPC

}

# eQTL saturation
p5 <- {
  qtl.rate %>% 
    ggplot(., 
           aes(x=maxPC, 
               y=n)) +
    geom_line(color="blue", linewidth=0.8) +
    geom_point(data=qtl.rate%>% filter(maxPC %% 5 == 0), size=2) +
    geom_vline(xintercept = best_maxPC, linetype=2, color="red") +
    labs(x="Max RNA PCs", y="Number of Significant QTLs") +
    themei
}

# Rate of change
p6 <- {
  qtl.rate %>%
    ggplot(., 
           aes(x=maxPC, 
               y=roc)) +
    geom_line() +
    geom_vline(xintercept = best_maxPC, linetype=2, color="red") +
    geom_hline(yintercept = 0) +
    labs(x="Max RNA PCs", y="Change in Significant QTLs") +
    themei
}

# Saturation point
pp <- (p5 |p6)
ggsave(plot = pp, paste0(prefix, ".saturation.pdf"), width = 16, height = 8)
#ggsave(plot = pp, paste0("/standard/vol185/cphg_Manichaikul/users/csm6hg/working_hf/figures/", prefix, ".saturation.pdf"), width = 16, height = 8)

# Output best k
write.table(
  best_maxPC,
  file = paste0("best_k_", prefix),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  sep = '\t')
