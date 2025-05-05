# Connor Murray
# Started 12.2.2024; modifed 2.10.2025
# analyzing TOPchef eQTLs colocalization
# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

# Libraries
library(data.table)
library(tidyverse)
library(argparse)

# Argument parser
parser <- ArgumentParser()
parser$add_argument("--wd", required=TRUE, help="The coloc output files directory.")
args <- parser$parse_args()

### Datasets & Setup ###

#wd="/standard/vol185/cphg_Manichaikul/users/csm6hg/nextflow_dna/output/coloc/"

# Coloc files
coloc <- list.files(path = args$wd, pattern = ".txt", full.names = T)
# coloc <- list.files(path = wd, pattern = ".txt", full.names = T)
coloc <- coloc[!coloc%like%"chrX"]
coloc <- coloc[coloc%like%"_gwas_HF_"]

# Read in streamlined GTF
gtf <- data.table(read_rds("/standard/vol185/cphg_Manichaikul/users/csm6hg/genome_files/gencode.v34.GRCh38.ERCC.genes.collapsed.streamlined.RDS"))

# Coloc run 
dt1 <- rbindlist(lapply(coloc, function(t) {fread(t, header = T)}), fill = T) %>% 
  left_join(gtf %>% 
              dplyr::select(-c(file, V7, V9)), 
            by=c("gene"="gene_edit","chrom")) %>% 
  mutate(dir=case_when(beta_qtl > 0 & beta_gwas > 0 ~ "+",
                       beta_qtl < 0 & beta_gwas < 0 ~ "-",
                       TRUE ~ "Flip"))

# Candidate colocalized genes !
candy <- dt1

# Define the significance threshold
threshold <- 0.8

# Define the custom theme
themei <- theme_bw() + 
  theme(axis.title.x = element_text(face = "bold", size = 18),
        axis.text.x = element_text(face = "bold", size = 18),
        axis.title.y = element_text(face = "bold", size = 18),
        axis.text.y = element_text(face = "bold", size = 18),
        legend.text = element_text(face = "bold", size = 18),
        legend.title = element_text(face = "bold", size = 18))

# Build plotting DF
p <- candy %>%
  filter(PP.H4 >= threshold) %>%
  dplyr::select(common_gene, gwas_pre, PP.H4, beta_qtl, dir) %>%
  distinct()

# Precompute vertical segment endpoints
segment_df <- p %>%
  group_by(common_gene) %>%
  summarise(min_y = min(gwas_pre), max_y = max(gwas_pre), .groups="drop")

# Make summary plot
plot1 <- {
  
  p %>%  
    ggplot(., 
           aes(x = common_gene, 
               y = gwas_pre)) +
    geom_segment(data = segment_df,
                 aes(x = common_gene, 
                     xend = common_gene, 
                     y = min_y, 
                     yend = max_y),
                 linetype = "dashed", color = "gray50", linewidth = 1) +
    geom_point(aes(shape = dir, fill= abs(beta_qtl)), size=4) +
    coord_flip() +
    scale_fill_viridis_c(name = "eQTL Effect Size") + 
    scale_shape_manual(name = "GWAS/eQTL Direction", 
                       values = c("+" = 24, "-" = 25, "Flip" = 21)) +
    themei +
    labs(x = "Colocalized Gene",
         size="eQTL Effect",
         y = "GWAS Dataset") +
    theme(legend.position = "right", 
          axis.text.x = element_text(face = "bold.italic", size = 14, angle = 45, hjust = 1),
          axis.text.y = element_text(face = "bold.italic", size = 12),
          axis.title = element_text(face = "bold", size = 14),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 14))

}

# Save output image
ggsave(plot = plot1, filename = "coloc_genes.pdf", dpi = 300, width = 12, height = 14)

# Output results
write_delim(candy, file = "coloc_eqtl_candidates_full.txt", delim = "\t")

# List to followup on with LD
candyi <- data.table(candy %>% select(chrom, minPos, maxPos, common_gene) %>% distinct())

write_delim(candyi, file = "coloc_eqtl_candidates_regions.txt", delim = "\t", col_names = F)
