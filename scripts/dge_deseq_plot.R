#!/usr/bin/env Rscript

# Differential Gene Expression Analysis with argparse and Enhanced Plotting
# By: Connor Murray (modified)
# Started 10.17.2024; modified 2.24.2025

# Libraries
library(argparse)
library(tidyverse)
library(data.table)
library(DESeq2)
library(apeglm)

# Argument Parsing
parser <- ArgumentParser(description='Differential gene expression analysis script with enhanced plotting')
parser$add_argument("--metadata", type="character", default="/standard/vol185/cphg_Manichaikul/users/csm6hg/metadata/metadata_10_17_2024_CSM.txt",
                    help="Metadata file")
parser$add_argument("--gtf", type="character", default="/standard/vol185/cphg_Manichaikul/users/csm6hg/genome_files/gencode.v34.GRCh38.ERCC.genes.collapsed.streamlined.RDS",
                    help="Genome GTF file streamlined")
parser$add_argument("--dge_list", type="character", default="/standard/vol185/cphg_Manichaikul/users/csm6hg/nextflow_dna/output/dge/sig_dge_DCM.csv",
                    help="DGE Genes")
args <- parser$parse_args()

# -------------------------------
# Input Files
# -------------------------------
meta <- data.table(fread(args$metadata, header = TRUE))
gtf <- data.table(readRDS(args$gtf))
res <- data.table(fread(args$dge_list, header = TRUE)) %>% 
  left_join(gtf %>% 
              select(gene, gene_edit, common_gene), 
            by=c("gene_id"="gene"))

# -------------------------------
# Prepare Metadata and Sample Selection
# -------------------------------
meta.simple <- data.table(meta %>%  
                  filter(diagnosis_simple %in% 
                           c("ICM", "IDCM", "Non-Failing")) %>% 
                  select(sample = SAMPLE_ID_TOR, 
                         condition = diagnosis_simple))

# Remove missing sample
meta.simple <- meta.simple[sample != "TOR238072"]

# Prepare data for volcano plot
res1 <- data.table(data.frame(res) %>%
                     mutate(pa = -log10(padj),
                            gene = rownames(res)) %>%
                     left_join(gtf, by = c("gene" = "gene")))

# Basic volcano plot with color for significance
vplot <- {
  res1 %>%
    mutate(DE = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Upregulated",
      padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
      TRUE ~ "NS")) %>%
    ggplot(aes(x = log2FoldChange, y = pa)) +
    geom_point(aes(color = DE), alpha = 0.5, size = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 1) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", linewidth = 1) +
    scale_color_manual(values = c("Upregulated" = "cyan4",
                                  "Downregulated" = "purple3",
                                  "NS" = "gray50")) +
    labs(x = expression("log2(Fold Change)"),
         y = expression("-log10(FDR " * italic(p) * "-value)"),
         color = "DCM Expression") +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
          axis.title.x = element_text(face = "bold", size = 18),
          axis.title.y = element_text(face = "bold", size = 18),
          axis.text.x = element_text(face = "bold", size = 16),
          axis.text.y = element_text(face = "bold", size = 16),
          legend.title = element_text(face = "bold", size = 16),
          legend.text = element_text(face = "bold", size = 14), 
          legend.position = "right",
          aspect.ratio = 1)
}

ggsave(filename = args$volcano_plot, plot = volcano_plot, width = 8, height = 6)

# Common ggplot theme for consistency
themei <- theme_bw() + 
  theme(axis.title.x = element_text(face = "bold", size = 16),
        axis.text.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 16),
        axis.text.y = element_text(face = "bold", size = 14),
        legend.text = element_text(face = "bold", size = 14),
        legend.title = element_text(face = "bold", size = 16))

# Merge gene metadata into a test dataset (for plotting) using the normalized data
dt.test <- data.table(dt[Name %in% filt_norm_counts$Name] %>% 
                        mutate(gene_edit = sub("\\..*", "", Name)) %>% 
                        left_join(gtf, by = c("gene_edit" = "gene_edit")))

# Calculate distribution statistics across the selected samples (assuming count columns are those matching sample IDs)
dt.test[, `:=`(
  gene_avg = -log10(rowMeans(.SD, na.rm = TRUE)),
  gene_avg_slen = -log10(rowMeans(.SD, na.rm = TRUE) / len),
  gene_sd = apply(.SD, 1, sd, na.rm = TRUE)
), .SDcols = patterns("TOR")]

# Standard deviation histogram plot
p1 <- dt.test %>% 
  ggplot(aes(x = gene_sd)) +
  geom_histogram(bins = 100, fill = "steelblue", color = "black") +
  labs(title = "Distribution of Gene Count Standard Deviations",
       x = "Standard Deviation",
       y = "Number of Genes") +
  themei
ggsave(filename = args$sd_histogram, plot = p1, width = 8, height = 6)

# Additional plots (e.g., gene length vs. average expression)
p2 <- dt.test %>% 
  ggplot(aes(x = -log10(len), y = gene_avg)) +
  geom_point(alpha = 0.4, color = "darkgreen") +
  geom_smooth(method = "lm", color = "red") +
  labs(title = "Gene Length vs. Average Expression",
       x = "-log10(Gene Length)",
       y = "-log10(Gene Count Average)") +
  themei
ggsave(filename = "gene_length_vs_expression.png", plot = p2, width = 8, height = 6)

# Linear model summary printed to console
lm_fit <- lm(-log10(len) ~ gene_avg, data = dt.test)
print(summary(lm_fit))

