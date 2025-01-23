# Connor Murray
# Started 10.21.2024
# analyzing TOPchef eQTLs 

# Libraries
library(tidyverse)
library(stringi)
library(data.table)
library(patchwork)
library(future.apply)

# Set up parallel backend, using all available cores
plan(multisession, workers = 4)

# Working directory 
setwd("/standard/vol185/cphg_Manichaikul/users/csm6hg")

# Metadata
meta <- data.table(fread("metadata/metadata_10_17_2024_CSM.txt", header = T))

# Candidate genes - DCM
dcm_cand <- data.table(fread("data/Gene-list-for-eQTL-analysis-SNC-10-9-2024.csv",
                             header = T))

# eQTL results
qtl.files <- list.files("/standard/vol185/cphg_Manichaikul/users/csm6hg/nextflow_dna/output/tensorqtl", 
                        pattern = "cis_qtl.txt.gz", full.names = T)

# List reading files
qtl <- rbindlist(future_lapply(qtl.files, function(file) {
  #data <- arrow::read_parquet(file)
  data <- fread(file)
  data[, .id := basename(file)]  # Add the file name as the .id column
  #print(basename(file))
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
gene.gtf <- data.table(readRDS("genome_files/gencode.v34.GRCh38.ERCC.genes.collapsed.streamlined.RDS"))

# Add gene metadata
qtl.dt <- data.table(qtl.dt %>% 
             left_join(gene.gtf %>% 
                      select(-c(file,V9)), 
                      by = c("phenotype_id" = "gene_edit", 
                             "chrom"="chrom")))


#qtl.dt[padj> -log10(0.05/21119)]
# saveRDS(qtl.dt, file = "data/qtl.rna.saturation.nomaf.11.12.24.rds")
#qtl.dt <- data.table(readRDS("data/qtl.rnapc.saturation.nomaf.11.12.24.rds"))

#### Analysis ####

# See candidate genes
qtl.dt[common_gene %in% dcm_cand$Gene][maxPC==11][pval_perm < cutoff]

# General stats
hist(qtl.dt$padj)
hist(qtl.dt$pval_nominal)
hist(qtl.dt$pval_perm)
hist(qtl.dt$slope)

# Common theme element
themei <- {
  theme_bw() + 
    theme(axis.title.x = element_text(family = "bold", size = 16),
          axis.text.x = element_text(family = "bold", size = 14),
          axis.title.y = element_text(family = "bold", size = 16),
          axis.text.y = element_text(family = "bold", size = 14),
          legend.text = element_text(family = "bold", size = 14),
          legend.title = element_text(family = "bold", size = 16)) 
}

# Pvalue cutoff
cutoff = 0.05

# qtl working
qtl.w <- data.table(qtl.dt)

# Distribution of p-values
p1 <- {
  
  qtl.w %>% 
    ggplot(., 
           aes(x=pval_perm)) +
    geom_histogram(bins=30) +
    geom_vline(xintercept = cutoff, 
               color="red") +
    labs(x="Permutation adj. p-value",
         y="Number of Genes") +
    theme_bw() +
    themei
  
}

p2 <- {
  
  qtl.w %>% 
    ggplot(., 
           aes(x=slope,
               y=-log10(pval_perm))) +
    geom_point(alpha=0.4) +
    geom_hline(yintercept = -log10(cutoff), 
               color="red") +
    labs(x="Effect Size",
         y="-log10(Perm. Adj. p-value)") +
    theme_bw() + 
    themei
  
}

p3 <- {
  
  qtl.w %>%
    mutate(sig=case_when(pval_perm < cutoff ~"Sig.",
                         TRUE ~ "Not Sig.")) %>% 
    ggplot(., 
           aes(x=af, 
               fill=sig)) +
    geom_histogram(position = "dodge") +
    labs(x="Allele Frequency",
         y="Number of Genes",
         fill="cis-eQTL") +
    theme_bw() + 
    themei
  
}

p4 <- {
  
  qtl.w %>%
    mutate(sig=case_when(pval_perm < cutoff ~"Sig.",
                         TRUE ~ "Not Sig.")) %>% 
    ggplot(., 
           aes(x=af,
               y=slope,
               color=sig)) +
    geom_point(alpha=0.4) + 
    labs(x="Allele Frequency",
         y="Effect Size",
         color="cis-eQTL") +
    theme_bw() + 
    themei
  
}

# Combine
pp <- ((p1 + p2) / (p3 + p4)) + plot_layout(guides = 'collect')

# Save
ggsave(plot = pp, 
       filename = "plots/qtl.sigPerm.maxPC11.png", 
       width = 10, height = 6)

### Saturation test

# Number of sig eQTLs and rate of change
qtl.rate <- data.table(qtl.dt %>% 
                 group_by(maxPC) %>% 
                 summarize(n = sum(pval_perm < cutoff)) %>% 
                 arrange(maxPC) %>%  
                 mutate(roc = (n - lag(n)) / (maxPC - lag(maxPC))))

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

# eQTL saturation
p5 <- {
  
  #qtl.dt[maxPC%in% c(2:32)] %>%
  qtl.rate %>% 
    ggplot(., 
           aes(x=maxPC, 
               y=n)) +
    geom_line(color="blue", size=0.8) +
    geom_point(data=qtl.rate%>% filter(maxPC %% 5 == 0), size=2) +
    geom_vline(xintercept = 11, linetype=2, color="red") +
    labs(x="Max RNA PCs", y="Number of Significant eQTLs") +
    themei 
  
}

# Rate of change
p6 <- {
  
  qtl.rate %>%
    ggplot(., 
           aes(x=maxPC, 
               y=roc)) +
    geom_line() +
    geom_vline(xintercept = 11, linetype=2, color="red") +
    geom_hline(yintercept = 0) +
    labs(x="Max RNA PCs", y="Change in Significant eQTLs") +
    themei 
  
}

# Extract top 100 genes for maxPC == 11
top_100_genes_maxPC11 <- qtl.dt %>%
  filter(maxPC == 11) %>%
  slice_max(pval_perm, n = 100) %>%
  pull(phenotype_id)

# Calculate the percentage of overlap for each maxPC (2 to 100) with maxPC == 11
op <- qtl.dt %>%
  group_by(maxPC) %>%
  slice_max(pval_perm, n = 100) %>%
  summarize(overlap = sum(phenotype_id %in% top_100_genes_maxPC11),
            percentage_overlap = (overlap / 100) * 100)

# Distribution of p-values
library(ggridges)
p7 <- {
  
  qtl.dt %>% 
    ggplot(.,aes(x=pval_perm, y=maxPC, group=maxPC)) +
    geom_density_ridges(data=qtl.dt %>% 
                          filter(maxPC %% 20 == 0) %>% 
                          left_join(qtl.rate), 
                        aes(x=pval_perm, 
                            y=maxPC, 
                            group=maxPC, 
                            fill=n), scale=1) +
    geom_vline(xintercept = cutoff, linetype=2, color="red") +
    viridis::scale_fill_viridis() +
    labs(x="Permutation adj. p-value)", 
         y="Max RNA PCs", 
         fill="Sig. eQTLs") +
    themei
  
}

# concordance of top 100 genes with PC11
p8 <- {
  
  op %>% 
    ggplot(., aes(x=maxPC, y=overlap)) +
    geom_line(size=1) +
    themei +
    labs(x="Max RNA PCs", y="Overlap of Top 100 eQTLs (%)", 
         title="Concordance of RNA PC 1-11") 
  
}

# correlation of significant genes
p9 <- {
  
  cor_list[!V1==V2][V1%in%c(1:40)][V2%in%c(1:40)] %>% 
    mutate(V1 = as.numeric(V1),
           V2 = as.numeric(V2)) %>% 
    ggplot(., aes(x=V1, y=V2, fill=N*N)) +
    geom_raster() +
    viridis::scale_fill_viridis(option = "B", begin = 1, end=0) +
    themei +
    labs(x="Max RNA PCs", 
         y="Max RNA PCs", 
         fill="R^2 of Sig. Genes") +
    coord_fixed()
  
}

# Step 1: Create cumulative positions
chrom_sizes <- qtl.dt %>%
  group_by(chrom) %>%
  summarise(chrom_max = max(stop))  # Get max position for each chromosome

# Add cumulative positions to qtl.dt
qtl.dt1 <- qtl.dt %>%
  left_join(chrom_sizes %>%
              mutate(cum_offset = cumsum(as.numeric(lag(chrom_max, default = 0)))),  # Cumulative offset for each chromosome
            by = "chrom") %>%
  mutate(cum_pos = cum_offset + (start + stop) / 2) %>%
  mutate(chrom_color = factor(as.numeric(str_remove(chrom, "chr")) %% 2)) %>% 
  mutate(chrom_color = factor(chrom_color, levels = as.character(1:22)))  

# gwas
p10 <- {
  
  qtl.dt1[maxPC==11] %>% 
    ggplot(., aes(x = cum_pos/1000000,
                  y = -log10(pval_nominal),
                  color = chrom)) +
    geom_point(alpha=0.4) +
    geom_hline(yintercept = -log10(0.05/21119), linetype=2, size=1) +
    themei +
    theme(legend.position = "none") +
    labs(x="Position (Mbp)", 
         y="-log10(eQTL p-value)", 
         color="")
  
}

# Saturation point
pp2 <- ((p5 / p6) | (p7 / p8)) / (p9 )

ggsave(plot = pp2, "eqtl.rnaconcord.png", width = 8, height = 8)















# Testing
qtl.dt[maxPC==11][pval_perm<0.0001]

# Some candidate genes
mean(qtl.dt[common_gene=="MFGE8"]$padj)
mean(qtl.dt[common_gene=="PSRC1"]$padj)
mean(qtl.dt[common_gene=="MYOZ1"]$padj)
