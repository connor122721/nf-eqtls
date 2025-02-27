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
library(pheatmap)

# Argument Parsing
parser <- ArgumentParser(description='Differential gene expression analysis script with enhanced plotting')
parser$add_argument("--metadata", type="character", default="/standard/vol185/cphg_Manichaikul/users/csm6hg/metadata/metadata_10_17_2024_CSM.txt",
                    help="Metadata file")
parser$add_argument("--gtf", type="character", default="/standard/vol185/cphg_Manichaikul/users/csm6hg/genome_files/gencode.v34.GRCh38.ERCC.genes.collapsed.streamlined.RDS",
                    help="Genome GTF file streamlined")
parser$add_argument("--dge_list", type="character", default="/standard/vol185/cphg_Manichaikul/users/csm6hg/nextflow_dna/output/dge/sig_dge_DCM.csv",
                    help="DGE Genes")
parser$add_argument("--counts", type="character", default="/standard/vol185/cphg_Manichaikul/users/csm6hg/nextflow_dna/output/rna/norm_medrat.tsv",
                    help="Normalized gene expression")
parser$add_argument("--outdir", type="character", default="/standard/vol185/cphg_Manichaikul/users/csm6hg/nextflow_dna/output/dge",
                    help="Output directory")
args <- parser$parse_args()

# Input Files
meta <- data.table(fread(args$metadata, header = TRUE))
gtf <- data.table(readRDS(args$gtf))
res <- data.table(fread(args$dge_list, header = TRUE)) %>%
  select(-c(common_gene, gene)) %>% 
  mutate(gene=tstrsplit(gene_id, ".", fixed=T)[[1]]) %>% 
  left_join(gtf %>% select(gene=gene_edit, common_gene), 
            by=c("gene"))
outdir=args$outdir
mat <- data.table(fread(args$counts)) %>% 
  left_join(meta, by=c("V1"="SAMPLE_ID_TOR"))
gene_cols <- setdiff(colnames(mat), "V1")
gene_ids <- tstrsplit(gene_cols, ".", fixed = TRUE)[[1]]
keep_idx <- which(gene_ids %in% gtf$gene_edit)
keep_cols <- c("V1", gene_cols[keep_idx])
mat <- mat[, ..keep_cols]
new_gene_names <- gtf$common_gene[match(gene_ids[keep_idx], gtf$gene_edit)]

# Rename the gene columns (the first column remains unchanged)
setnames(mat, old = gene_cols[keep_idx], new = new_gene_names)

# Prepare Metadata and Sample Selection
meta.simple <- data.table(meta %>%  
                  filter(diagnosis_simple %in% 
                           c("ICM", "IDCM", "Non-Failing")) %>% 
                  select(sample = SAMPLE_ID_TOR, 
                         condition = diagnosis_simple,
                         sex=Gender))

# Remove missing sample
meta.simple <- meta.simple[sample != "TOR238072"]
meta.simple <- meta.simple[!sex==""]

# Prepare data for volcano plot
res1 <- data.table(data.frame(res) %>%
                     mutate(pa = -log10(padj),
                            gene = rownames(res)))

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
    labs(x = expression(bold("log2(Fold Change)")),
         y = expression(bold("-log10(FDR " * italic(p) * "-value)")),
         color = "DCM Expression:") +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
          axis.title.x = element_text(face = "bold", size = 18),
          axis.title.y = element_text(face = "bold", size = 18),
          axis.text.x = element_text(face = "bold", size = 16),
          axis.text.y = element_text(face = "bold", size = 16),
          legend.title = element_text(face = "bold", size = 16),
          legend.text = element_text(face = "bold", size = 14), 
          legend.position = "top",
          aspect.ratio = 1)
  
}

ggsave(path = outdir, 
       filename = "DCM_volcano.pdf", 
       plot = vplot, width = 8, height = 8)

# Make heatmap
res.dt <- data.table(res %>%
            mutate(DE = case_when(padj < 0.05 & log2FoldChange > 1  ~ "Upregulated",
                                  padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
                                  TRUE ~ "NS")))

# Function to create heatmap of expression profiles
heatmap_dge <- function(genes, out_name) {

  cols <- intersect(colnames(mat), genes)
  samples <- meta.simple[meta.simple$condition %in% c("IDCM", "Non-Failing"), "sample"][[1]]
  matt <- subset(mat, V1 %in% samples, select = c("V1", cols))
  
  # Set sample names as row names and remove the V1 column
  rownames(matt) <- matt$V1
  matt$V1 <- NULL
  matt_t <- t(matt)
  colnames(matt_t) <- rownames(matt)
  
  # Prepare sample metadata; row names should match the column names of the transposed matrix
  meta.matt <- meta.simple[meta.simple$sample %in% samples, ]
  metamm <- meta.matt %>% select("Case" = condition, Sex=sex)
  rownames(metamm) <- meta.matt$sample
  
  # Optionally, define custom colors for the condition factor
  condition_colors <- c("IDCM" = "cyan3", "Non-Failing" = "pink3")
  sex_colors <- c("Female" = "lightgreen", "Male" = "green4")
  ann_colors <- list(Case = condition_colors, Sex = sex_colors)
  
  # Run pheatmap with the transposed matrix and sample annotations for condition only
  pheatmap(matt_t, 
           cluster_rows = T, 
           show_colnames = F, 
           annotation_names_col = F,
           show_rownames = T,
           annotation_col = metamm,
           annotation_colors = ann_colors,
           fontsize = 14, 
           scale = "row",
           fontsize_row = 14, 
           height = 25,
           width = 20,
           filename = paste0(outdir,out_name))

}

# Subset matrix to top upregulated genes
up_genes <- res.dt[DE == "Upregulated" & log2FoldChange>2.2, common_gene]
down_genes <- res.dt[DE == "Downregulated" & log2FoldChange< -2.2, common_gene]

heatmap_dge(up_genes,
            "DCM_upreg_heat.pdf")

heatmap_dge(down_genes,
            "DCM_downreg_heat.pdf")

####################################
## Gene orthology enrichment for upexpression ##
####################################

# Enrichment analysis function
library(clusterProfiler)
library(org.Hs.eg.db)
library(foreach)
library(DOSE)

# Universe: the set of longest transcript proteins (unique)
univ <- unique(res.dt$common_gene)

# Clusterprofiler function
go.test.clusterpro <- function(cand, out_name) {
  #cand=data.table(res.dt %>% filter(DE=="Upregulated" & log2FoldChange > cut) %>% dplyr::select(common_gene))$common_gene; out_name="DCM_upregulated"
  
  # Enrichment test using clusterProfiler's enricher()
  go.test <- enrichGO(gene = cand,
                      OrgDb = org.Hs.eg.db,
                      pvalueCutoff = 0.05,
                      keyType = "SYMBOL",
                      pAdjustMethod = "fdr",
                      minGSSize = 10,
                      maxGSSize = 500,
                      qvalueCutoff = 0.01,
                      universe = univ)

  go.test2 <- data.table(go.test@result) %>% filter(p.adjust < 0.05)
  
  # Iterate over each enriched GO term to compute odds ratios and add term descriptions
  pp <- foreach(i = seq_len(nrow(go.test2)), .combine = "rbind", .errorhandling = "remove") %do% {
    current_id <- go.test2$ID[i]
    message(paste(i, current_id, sep = " | "))
    go.test2[ID == current_id] %>% 
      mutate(term = go2term(ID)$Term,          # Make sure go2term() is defined
             x = eval(parse(text = GeneRatio)),
             y = eval(parse(text = BgRatio)),
             odd = x / y)
  }
  
  # Plot enriched GO terms
  plot.go <- { 
    
    pp %>% 
      filter(p.adjust < 0.05) %>% 
      arrange(desc(p.adjust)) %>% 
      ggplot(aes(x = log2(odd),
                 y = reorder(term, log2(odd)),
                 size = Count,
                 color = log10(p.adjust))) +
      geom_point() +
      labs(x = expression(bold("log2(Enrichment)")),
           y = "",
           color = expression(bold("-log10(FDR " * italic(p) * ")")),
           size = "Candidate gene number") +
      theme_bw() +
      viridis::scale_color_viridis(option = "magma") +
      theme(title = element_text(face = "bold", size = 15),
            legend.text = element_text(face = "bold", size = 14),
            legend.title = element_text(face = "bold", size = 16),
            legend.background = element_blank(),
            axis.text.x = element_text(face = "bold", size = 15),
            axis.text.y = element_text(face = "bold", size = 15),
            axis.title.x = element_text(face = "bold", size = 18),
            axis.title.y = element_text(face = "bold", size = 18))
  }
  
  ggsave(plot = plot.go, 
         path = outdir,
         filename = paste0(out_name, ".goterm.pdf"), 
         width = 20, height = 12)
  
  return(go.test)

}

# Effect size cutoff
cut=1.5

# Upregulated!
up <- go.test.clusterpro(data.table(res.dt %>% 
                     filter(DE=="Upregulated" & log2FoldChange > cut) %>% 
                     dplyr::select(common_gene))$common_gene, 
                     "DCM_upregulated")

# Downregulated!
down <- go.test.clusterpro(data.table(res.dt %>% 
                     filter(DE=="Downregulated" & log2FoldChange < -cut) %>% 
                     dplyr::select(common_gene))$common_gene, 
                     "DCM_downregulated")

######################
## Pathway analysis ##
######################
library(enrichplot)
library(patchwork)

# Create a gene-concept network
upnet <- cnetplot(up, 
         color_gene='steelblue', 
         showCategory = 5) +
  labs(title="DCM Upregulated")

downet <- cnetplot(down, 
         color_gene='steelblue', 
         showCategory = 5) +
  labs(title="DCM Downregulated")

fin <- (downet | upnet) 

ggsave(path = outdir, 
       filename = "DCM_network.pdf", 
       plot = fin, width = 12, height = 10)
