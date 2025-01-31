# Connor Murray
# Started 12.2.2024
# analyzing TOPchef eQTLs colocalization
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
parser$add_argument("--wd", required=TRUE, help="The coloc output files directory.")
args <- parser$parse_args()

### Datasets & Setup ###

# Coloc files
coloc <- list.files(path = args$wd, pattern = ".txt", full.names = T)
coloc <- coloc[!coloc%like%"chrX"]
coloc <- coloc[coloc%like%"_gwas_HF_"]

# Read in streamlined GTF
gtf <- data.table(read_rds("/standard/vol185/cphg_Manichaikul/users/csm6hg/genome_files/gencode.v34.GRCh38.ERCC.genes.collapsed.streamlined.RDS"))

# Coloc run 
dt1 <- rbindlist(lapply(coloc, function(t) {fread(t, header = T)})) %>% 
  left_join(gtf %>% 
              select(-c(file, V7, V9)), 
            by=c("gene"="gene_edit","chrom"))

# Candidate colocalized genes !
candy <- dt1[PP.H4 >= 0.5]

### Plot ###

# Transform the data to wide format with PP.H4 scores
dt.co <- data.table(candy %>% 
                      select(common_gene, gwas, PP.H4) %>% 
                      pivot_wider(values_from = PP.H4, 
                                  names_from = gwas))

# Define the significance threshold
threshold <- 0.8

# Identify significant genes for each GWAS and retain PP.H4 scores or set to NA
dt.co[, `:=`(
  Levin2022 = ifelse(levin22_gwas_HF > threshold, levin22_gwas_HF, NA),
  Shah2020 = ifelse(shah20_gwas_HF > threshold, shah20_gwas_HF, NA))]

# Pivot the data to long format for plotting
b <- data.frame(dt.co %>% 
                  select(common_gene, Levin2022, Shah2020) %>%
                  pivot_longer(cols = c(Levin2022, Shah2020),
                               names_to = "GWAS",
                               values_to = "PP.H4"))

# Define the custom theme
themei <- theme_bw() + 
  theme(axis.title.x = element_text(face = "bold", size = 18),
        axis.text.x = element_text(face = "bold", size = 18),
        axis.title.y = element_text(face = "bold", size = 18),
        axis.text.y = element_text(face = "bold", size = 18),
        legend.text = element_text(face = "bold", size = 18),
        legend.title = element_text(face = "bold", size = 18))

# Create the Upset-like plot with point sizes based on PP.H4 scores
plot1 <- {
  b %>% 
    filter(!is.na(PP.H4)) %>% 
    group_by(common_gene) %>% 
    mutate(min_y = min(GWAS, na.rm = TRUE),
           max_y = max(GWAS, na.rm = TRUE)) %>% 
    ggplot(aes(x = common_gene, y = GWAS)) +
    geom_segment(aes(x = common_gene, 
                     xend = common_gene, 
                     y = min_y, 
                     yend = max_y), 
                 linetype = "dashed", 
                 color = "gray", 
                 linewidth = 2) +
    geom_point(aes(color = PP.H4), size=8) +
    scale_color_continuous(type = "viridis") +
    themei +
    labs(x = "Colocalized Gene",
         y = "Query Dataset", 
         title = "Colocalization with TOPCHef eQTLs (PP.H4 > 0.8)") +
    theme(legend.position = "right",
          axis.text.x = element_text(face = "bold.italic", size = 18, angle = 45, hjust = 1),
          axis.text.y = element_text(face = "bold.italic", size = 18),
          axis.title = element_text(face = "bold", size = 18),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 18))
}

# Save output image
ggsave(plot = plot1, filename = "coloc_genes.pdf", dpi = 300, width = 13, height = 7)

# Output results
write_delim(candy, file = "coloc_eqtl_candidates_full.txt", delim = "\t")

# List to followup on with LD
candyi <- data.table(candy %>% select(chrom, minPos, maxPos, common_gene) %>% distinct())

write_delim(candyi, file = "coloc_eqtl_candidates_regions.txt", delim = "\t", col_names = F)
