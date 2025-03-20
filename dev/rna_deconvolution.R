# Connor S. Murray
# Bulk RNA-seq deconvolution for heart tissue
# Started: 3.11.2025

# Libraries
library(MuSiC)
library(Seurat)
library(data.table)
library(tidyverse)
library(zellkonverter)
library(SeuratDisk)
library(Biobase)
library(SummarizedExperiment)

# Read in scRNAseq dataset from heart tissue
d = "/standard/vol185/cphg_Manichaikul/users/csm6hg/data/Global_raw.h5ad"
# d = "/standard/vol185/cphg_Manichaikul/users/csm6hg/data/visium-OCT_LV_lognormalised.h5ad"
sc <- zellkonverter::readH5AD(d) # Reads in SingleCellExperiment
names(sc@assays) <- "counts"
set.seed(100)

# Convert to Seurat input
#sce <- as.Seurat(sc, counts = "X", data = NULL)
# Normalize, find variable features, scale data, and run PCA
#sce <- NormalizeData(sce)
#sce <- FindVariableFeatures(sce)
#sce <- ScaleData(sce)
#sce <- RunPCA(sce)
#sce <- RunUMAP(sce, dims = 1:10)
#DimPlot(sce, reduction = "umap", label = TRUE) + ggtitle("UMAP before filtering")

# Pt 2. Prepare Bulk RNA-seq

# Read in GTF
gtf <- data.table(readRDS("/standard/vol185/cphg_Manichaikul/users/csm6hg/genome_files/gencode.v34.GRCh38.ERCC.genes.collapsed.streamlined.RDS"))

# Read in normalized TMM RNAseq data
norm <- "/standard/vol185/cphg_Manichaikul/users/csm6hg/nextflow_dna/output/rna/norm_tmm.tsv"
norm_tmm <- data.table(fread(norm, header = T)) %>% 
  mutate(gene = tstrsplit(Name, ".", fixed=T)[[1]]) %>% 
  left_join(gtf %>% select(gene=gene_edit, common_gene))

raw <-"/standard/vol185/TOPMed/TOPCHef/82214/topmed-dcc/exchange/phs002038_TOPMed_TOPCHeF/Omics/RNASeq/release3/TOPMed_Taylor_P4.RNASeQCv2.3.6_gene_reads.gct.gz"
raw <- data.table(fread(raw, header = T))
raw <- raw[!duplicated(Description)] %>% 
  filter(Name %in% norm_tmm$Name) %>% 
  mutate(gene=tstrsplit(Name, ".", fixed=T)[[1]])
# raw <- raw[1:100] # TESTING
m <- as.matrix(raw %>% select(-c(Name, gene, Description)))
row.names(m) <- raw$gene
m <- ExpressionSet(assayData = m)
bulk.mtx = exprs(m)

# Change rownames to gene symbols
norm_tmm <- norm_tmm[!duplicated(common_gene)]
# norm_tmm <- norm_tmm[1:100] # TESTING

m_norm <- as.matrix(norm_tmm %>% select(-c(Name, gene, common_gene)))
#row.names(m) <- norm_tmm$common_gene
row.names(m_norm) <- norm_tmm$gene
m_norm <- ExpressionSet(assayData = m_norm)
bulkNorm.mtx = exprs(m_norm)

# Get the common gene names between the single-cell data and the bulk data
common_genes <- intersect(rownames(sc), rownames(m))
cat("Number of common genes:", length(common_genes), "\n")

# Subset the single-cell dataset to keep only those common genes
sc_sub <- sc[common_genes, ]

# Adjust the threshold for n_genes as needed.
cells_to_keep <- (colData(sc_sub)$n_genes >= 1000 & 
                  colData(sc_sub)$total_counts >= 800 &
                  colData(sc_sub)$modality == "scRNA")
sc_filtered <- sc_sub[, cells_to_keep]
sample_meta <- unique(data.frame(sangerID = colData(sc_filtered)$sangerID,
                                 gender = colData(sc_filtered)$gender))

# Randomly select 5 samples per gender
selected_samples <- sample_meta %>%
  group_by(gender) %>%
  sample_n(5) %>%
  pull(sangerID)

cat("Selected samples:", paste(selected_samples, collapse = ", "), "\n")
sc_subset <- sc_filtered[, colData(sc_filtered)$sangerID %in% selected_samples]

# Pt 3. Run Deconvolution with MuSiC
deconv_results <- music_prop(bulk.mtx = bulk.mtx,
                             sc.sce = sc_subset,
                             clusters = "cell_type",
                             samples = "sangerID",
                             verbose = TRUE, 
                             select.ct = c("Atrial Cardiomyocyte",
                                           "Fibroblast",
                                           "Mast cell", 
                                           "Neural cell",
                                           "Endothelial cell",
                                           "Myeloid",
                                           "Ventricular Cardiomyocyte"))

# Metadata
meta <- fread("/standard/vol185/cphg_Manichaikul/users/csm6hg/metadata/metadata_10_17_2024_CSM.txt")

# Convert the matrix to a data frame and add a column for cell types
prop_mat <- deconv_results$Est.prop.weighted
prop_df <- as.data.frame(prop_mat)
prop_df$sample <- rownames(prop_df)

# Convert to long format so that each row represents a sample-cell type combination
prop_long <- data.table(pivot_longer(prop_df, 
                                     cols = -sample, 
                                     names_to = "cell_type", 
                                     values_to = "proportion") %>% 
                right_join(meta %>%
                    filter(diagnosis_simple %in% c("IDCM", "ICM", "Non-Failing")) %>% 
                    select(sample=SAMPLE_ID_TOR, diagnosis_simple, Gender, Age_at_collection)))

# Combine verbose cell type names into simpler categories
prop_long <- prop_long %>%
  mutate(cell_type = case_when(
    #cell_type %in% c("Atrial Cardiomyocyte", "Ventricular Cardiomyocyte") ~ "Cardiomyocytes",
    cell_type == "Fibroblast" ~ "Fibroblast",
    cell_type == "Ventricular Cardiomyocyte" ~ "Cardiomyocyte",
    cell_type == "Atrial Cardiomyocyte" ~ "Cardiomyocyte",
    #cell_type %in% c("Mast cell", "Myeloid") ~ "Blood cells",
    #cell_type == "Neural cell" ~ "Neurons",
    TRUE ~ cell_type))

# Some summaries: average proportion by simplified cell type and diagnosis
prop_long %>% 
  group_by(cell_type, diagnosis_simple) %>% 
  summarize(m = mean(proportion, na.rm = TRUE))

# Define the custom theme
themei <- theme_bw() + 
  theme(axis.title.x = element_text(face = "bold", size = 18),
        axis.text.x = element_text(face = "bold", size = 18),
        axis.title.y = element_text(face = "bold", size = 18),
        axis.text.y = element_text(face = "bold", size = 18),
        legend.text = element_text(face = "bold", size = 18),
        legend.title = element_text(face = "bold", size = 18))

# Create the stacked bar plot
desired_order <- c("Cardiomyocyte", "Endothelial cell", 
                   "Fibroblast", "Myeloid", "Neural cell")
prop_long$cell_type <- factor(prop_long$cell_type, levels = desired_order)

p1 <- {
  prop_long[!Gender == ""] %>% 
    ggplot(aes(x = proportion*100, 
               y = reorder(sample, Age_at_collection),
               fill = cell_type)) +
    geom_bar(stat = "identity") +
    facet_wrap(~paste(diagnosis_simple,Gender), scales = "free_y", ncol = 2) +
    theme_bw() +
    labs(x = "Cell type proportions (%)",
         y = "Sample",
         fill = "") +
    themei +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          aspect.ratio = 1,
          strip.text = element_text(face = "bold", size = 18))
}

p2 <- {
  na.omit(prop_long) %>% 
    ggplot(., 
           aes(x = diagnosis_simple,
               y = proportion*100, 
               fill = diagnosis_simple)) +
    geom_violin() +
    geom_boxplot(fill="white", outlier.shape = NA) +
    scale_fill_manual(values = c("steelblue4", "blue2", "lightblue1")) +
    facet_wrap(~cell_type, 
               scales = "free_y", 
               ncol = 1) +
    theme_classic() +
    labs(x = "",
         y = "Cell type proportions (%)",
         fill="") +
    themei +
    theme(legend.position = "none",
          aspect.ratio = 1,
          axis.text.x = element_text(angle = -30),
          strip.text = element_text(face = "bold", size=12))
}

# Combine and output
library(patchwork)
p3 <- (p1+p2)

ggsave("/standard/vol185/cphg_Manichaikul/users/csm6hg/nextflow_dna/output/deconvolution_new.pdf", 
       p3, width = 16, height = 16)

saveRDS(prop_long, file = "/standard/vol185/cphg_Manichaikul/users/csm6hg/nextflow_dna/output/deconvolution.rds")
