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
parser$add_argument("--metadata", type="character", default="metadata/metadata_10_17_2024_CSM.txt",
                    help="Metadata file")
parser$add_argument("--gtf", type="character", default="genome_files/gencode.v34.GRCh38.ERCC.genes.collapsed.streamlined.RDS",
                    help="Genome GTF file streamlined")
parser$add_argument("--counts", type="character", default="/standard/vol185/TOPMed/TOPCHef/82214/topmed-dcc/exchange/phs002038_TOPMed_TOPCHeF/Omics/RNASeq/release3/TOPMed_Taylor_P4.RNASeQCv2.3.6_gene_reads.gct.gz",
                    help="Raw gene counts file")
parser$add_argument("--fingenes", type="character", default="qtl/run_10_8_24/filt_rnaseq_norm_topchef.bed",
                    help="Final genes used for eQTL")
parser$add_argument("--output", type="character", default="norm_counts_filtered_csm.rds",
                    help="Output file for filtered normalized counts")
args <- parser$parse_args()


# -------------------------------
# Input Files
# -------------------------------
meta <- data.table(fread(args$metadata, header = TRUE))
gtf <- data.table(readRDS(args$gtf))
fingenes <- data.table(fread(args$fingenes))
fingenes[, gene_edit := sub("\\..*", "", gene)]

# -------------------------------
# RNAseq Data & Raw Counts
# -------------------------------
dt <- data.table(fread(file.path(args$counts)))
cols <- colnames(dt)[c(1, which(colnames(dt) %in% meta$SAMPLE_ID_TOR))]

# -------------------------------
# Prepare Metadata and Sample Selection
# -------------------------------
meta.simple <- data.table(meta %>%  
                  filter(diagnosis_simple %in% c("ICM", "IDCM", "Non-Failing")) %>% 
                  select(sample = SAMPLE_ID_TOR, 
                         condition = diagnosis_simple))

# Remove missing sample
meta.simple <- meta.simple[sample != "TOR238072"]

# Randomly sample 20 individuals per condition
simp <- data.table(meta.simple %>% 
                     group_by(condition) %>% 
                     sample_n(20))$sample

# -------------------------------
# Filter and Format Count Data
# -------------------------------
# Select only desired samples
dt1 <- dt %>% select(Name, all_of(simp))
row.names(dt1) <- dt1$Name

# Add gene_edit column for merging
dt1[, gene_edit := sub("\\..*", "", Name)]

# ---------------------------------
# Differential Expression Analysis
# ---------------------------------
dds <- DESeqDataSetFromMatrix(countData = dt1 %>% 
                                select(-c(Name, gene_edit)),
                              colData = meta.simple[sample %in% simp],
                              design = ~ condition)
dds <- DESeq(dds)
rownames(dds) <- dt1$Name
res <- results(dds)
print(summary(res, alpha = 0.05))

# Save dispersion plot as PNG (base plot)
png("disp_plot.png", width = 800, height = 600)
plotDispEsts(dds, main = "Dispersion plot")
dev.off()

# -------------------------------
# MA Plot and Volcano Plot Enhancements
# -------------------------------
# Standard MA plot
png("ma_plot.png", width = 800, height = 600)
plotMA(res, main = "MA Plot", ylim = c(-2, 2))
dev.off()

# Apply shrinkage for log2 fold changes and replot MA
resLFC <- lfcShrink(dds = dds, coef = 2, type = "apeglm")
png("ma_plot_shrink.png", width = 800, height = 600)
plotMA(resLFC, main = "MA Plot (shrunken LFC)", ylim = c(-2, 2))
dev.off()

# Prepare data for volcano plot
res1 <- data.table(data.frame(res) %>%
                     mutate(pa = -log10(padj),
                            gene = rownames(res)) %>%
                     left_join(gtf, by = c("gene" = "gene")))

# Basic volcano plot with color for significance
volcano_plot <- res1 %>% 
  ggplot(aes(x = log2FoldChange, y = pa)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_point(aes(color = padj < 0.05), alpha = 0.7) +
  scale_color_manual(values = c("black", "red"), name = "padj < 0.05") +
  labs(title = "Volcano Plot",
       x = "log2(Fold Change)",
       y = "-log10(Adjusted p-value)") +
  theme_minimal(base_size = 14)

res1
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

