# Connor Murray
# Started 12.10.2024; modifed 1.30.2025
# analyzing TOPchef eQTLs and perform colocalization
# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

# Libraries
library(data.table)
library(coloc)
library(tidyverse)
library(locuscomparer)
library(foreach)
library(argparse)

# Argument parser
parser <- ArgumentParser()
parser$add_argument("--gwas", required=TRUE, help="Processed GWAS summary statistics.")
parser$add_argument("--eqtl", required=TRUE, help="cis-eQTL statistics.")
parser$add_argument("--chromosome", required=TRUE, help="Current chromosome.")
parser$add_argument("--gene", required=TRUE, help="Current gene.")

args <- parser$parse_args()

gwas_input <- args$gwas
eqtl_inpt <- args$eqtl
chromosome <- args$chromosome
gene <- args$gene

### Datasets & Setup ###
gtf <- data.table(readRDS("/standard/vol185/cphg_Manichaikul/users/csm6hg/genome_files/gencode.v34.GRCh38.ERCC.genes.collapsed.streamlined.RDS")) %>% 
  select(-c(V7:V9, file), gene_length=len)
files <- list.files(path = "coloc_gene/output", pattern = ".txt", full.names = T)
rsid <- fread("data/rsid_output.csv") %>% select(rsid, variant_id)

# Common theme element
themei <- {
  theme_bw() + 
    theme(axis.title.x = element_text(face = "bold", size = 16),
          axis.text.x = element_text(face = "bold", size = 14),
          axis.title.y = element_text(face = "bold", size = 16),
          axis.text.y = element_text(face = "bold", size = 14),
          legend.text = element_text(face = "bold", size = 14),
          legend.title = element_text(face = "bold", size = 16)) 
}

#### MAKE LD aware plots ####

# Read in topchef eQTLs and restrict to colocalized genes
qtl <- data.table(arrow::read_parquet(eqtl)) %>% 
  filter(phenotype_id == gene)

# Read in GWAS Shah 2020
gwas <- data.table(readRDS(gwas))

# Merge
dt <- data.table(na.omit(qtl %>%
        left_join(gwas, 
              by=c("variant_id"="snpID_hg38"))) %>% 
        left_join(gtf, by=c("phenotype_id"="gene")) %>% # Add gene info
        left_join(dt1 %>% select("phenotype_id"=gene, maxSNP)))

# Read in LD matrix MYOZ1
ld.file <- list.files(path = "linkage", 
                      pattern = gene, 
                      full.names = T)

# Read and process LD file
readLD_process <- function(ldFile) {
  
  # Read in LD file
  ld <- data.table(fread(input = ldFile, header = T) %>% 
                     mutate(SNPA=paste("chr", CHR_A, ":", BP_A, "[b37]", sep=""),
                            SNPB=paste("chr", CHR_B, ":", BP_B, "[b37]", sep="")) %>% 
                     filter(SNPA %in% dt$variant_id & SNPB %in% dt$variant_id) %>% 
                     mutate(R = sqrt(R2)) %>% 
                     select(SNPA, SNPB, R, R2))
  
  # Get all unique SNPs
  snp_list <- unique(c(ld$SNPA, ld$SNPB))
  snp_list <- sort(snp_list)  # Optional sorting
  
  # Create an empty matrix
  n_snps <- length(snp_list)
  ld_matrix <- matrix(0, nrow = n_snps, ncol = n_snps,
                      dimnames = list(snp_list, snp_list))
  
  # Create a mapping of SNP names to matrix indices for fast assignment
  snp_indices <- setNames(seq_along(snp_list), snp_list)
  
  # Assign R values to the matrix using indices
  ld_matrix[cbind(snp_indices[ld$SNPA], snp_indices[ld$SNPB])] <- ld$R2
  ld_matrix[cbind(snp_indices[ld$SNPB], snp_indices[ld$SNPA])] <- ld$R2
  
  # Set diagonal values to 1
  diag(ld_matrix) <- 1
  
  # Convert the LD matrix to a data.table
  ld_dt <- as.data.table(as.table(ld_matrix))
  setnames(ld_dt, c("V1", "V2", "N"), c("SNPA", "SNPB", "R2"))

  # Finish
  return(ld_dt)
}

# Process
ld_dt <- data.table(readLD_process(ld.file))

# Merge with dt MYOZ1
myo <- na.omit(dt[common_gene=="MYOZ1"] %>% 
  left_join(ld_dt %>% 
              filter(SNPA=="chr10:73646383[b37]"), 
            by=c("variant_id"="SNPB")))

# Plot gwas and eqtl statistics
library(ggrepel)

# MYOZ1 
plot2 <- {
  myo %>% 
    ggplot(aes(x = -log10(pval_nominal),
               y = -log(p),
               color = R2)) +
    geom_point(size=2) +
    geom_point(data = myo[variant_id == "chr10:73646383[b37]"], 
               size = 5, shape = 8, color = "black") +
    geom_label_repel(data = myo[variant_id == "chr10:73646383[b37]"], 
                     aes(label = SNP), alpha = 0.8) +
    facet_wrap(~common_gene, scales = "free") +
    scale_color_gradient2(low = "blue", mid = "green", high = "red", 
                          midpoint = 0.5, name = expression(italic(R)^2)) +
    labs(x = expression("-log(cis-eQTL " * italic(p) * "-value)"),
         y = expression("-log(Shah 2020 GWAS " * italic(p) * "-value)")) +
    themei +
    theme(strip.text = element_text(face = "bold", size = 16),
          legend.title = element_text(face = "bold", size = 20))
}

ggsave(plot = plot1, 
       filename = "coloc_gene/myoz1_topchef_shah.coloc.png", 
       width = 8, height = 7)

plot2 <- {
  dt %>% 
    ggplot(., 
           aes(x=-log10(pval_nominal),
               y=-log(p),
               color=phenotype_id)) +
    geom_point(alpha=0.6) +
    geom_point(data=dt[common_gene=="AC073389.2"][variant_id=="chr10:73657491[b37]"], size=5, shape=8) +
    geom_point(data=myo[variant_id=="chr10:73646383[b37]"], size=5, shape=8) +
    geom_point(data=dt[common_gene=="SYNPO2L"][variant_id=="chr10:73661890[b37]"], size=5, shape=8) +
    geom_label_repel(data=dt[common_gene=="AC073389.2"][variant_id=="chr10:73657491[b37]"], aes(label=SNP), alpha=0.6) +
    geom_label_repel(data=dt[common_gene=="MYOZ1"][variant_id=="chr10:73646383[b37]"], aes(label=SNP), alpha=0.6) +
    geom_label_repel(data=dt[common_gene=="SYNPO2L"][variant_id=="chr10:73661890[b37]"], aes(label=SNP), alpha=0.6) +
    facet_wrap(~common_gene, scales = "free") +
    labs(x="-log(cis-eQTL nominal p-value)",
         y="-log(GWAS p-value)") +
    themei +
    theme(legend.position = "none",
          strip.text = element_text(family = "bold", size=16))

}

