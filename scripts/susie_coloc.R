# Connor Murray
# Started 11.18.2024; modified 2.11.2025
# analyzing TOPchef eQTLs and perform fine-mapping
# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

# Libraries
library(data.table)
library(susieR)
library(tidyverse)
library(argparse)
library(foreach)
library(purrr)

parser <- ArgumentParser()
parser$add_argument("--chromosome", required=TRUE, help="Chromosome to process, e.g. 'chr6'")
parser$add_argument("--coloc", required=TRUE, help="Path to coloc file using coloc.abf")
parser$add_argument("--qtl", required=TRUE, help="Path to QTL parquet file of nominal p-values")
parser$add_argument("--ld_dir", required=TRUE, help="Directory containing LD files")
parser$add_argument("--N", type="integer", default=516, help="Sample size (N1)")
args <- parser$parse_args()

chr <- args$chromosome
coloc_file <- args$coloc
qtl_file <- args$qtl
ld_dir <- args$ld_dir
N1 <- args$N
out_file <- args$out

# setwd("/standard/vol185/cphg_Manichaikul/users/csm6hg/nextflow_dna/"); chr="chr10"; coloc_file="output/coloc_no_se/coloc_eqtl_jurgens24_gwas_HF_chr10.txt";qtl_file="output/tensorqtl_nominal/topchef_chr10_MaxPC49.cis_qtl_pairs.chr10.parquet";ld_dir="output/linkage/"; N1=516;

### Datasets & Setup ###
gtf <- data.table(readRDS("/standard/vol185/cphg_Manichaikul/users/csm6hg/genome_files/gencode.v34.GRCh38.ERCC.genes.collapsed.streamlined.RDS")) %>% 
  select(-c(V7:V9, file), gene_length=len)

coloc <- fread(coloc_file) %>% left_join(gtf, by=c("gene"="gene_edit", "chrom"))
coloc <- coloc[min_p.eqtl < 1e-3 & PP.H4 > 0.01]
snps <- coloc[chrom %in% chr]

qtl <- data.table(arrow::read_parquet(qtl_file)) %>% 
  left_join(gtf, by=c("phenotype_id"="gene_edit"))
qtl[, snp := as.numeric(tstrsplit(variant_id, ":")[[2]])]

# Function to run susie with LD matrix input
susieWindow <- function(dt1, focalGene, N1) {
  #dt1=qtl; focalGene="BAG3"

  print(focalGene)
  
  # Define the window boundaries
  Min_pos <- snps[common_gene == focalGene]$minPos
  Max_pos <- snps[common_gene == focalGene]$maxPos
  
  # Filter dt1 to the desired window and gene
  dt_window <- dt1[snp >= Min_pos & snp <= Max_pos & common_gene == focalGene]
  
  # Ensure the variant IDs are character and unique (in case dt1 had duplicates)
  dt_window[, variant_id := as.character(variant_id)]
  dt_window <- unique(dt_window, by = "variant_id")

  # List the files for the given chromosome and gene pattern
  ld_files <- list.files(path = ld_dir, pattern = paste(chr, "_", focalGene, sep = ""), full.names = TRUE)
  ld_file <- ld_files[ld_files %like% "ld.gz"]
  snp_names <- fread(ld_files[ld_files %like% "snp_names.txt"], header = FALSE)
  
  # Read in the LD matrix (fread can handle gzipped files)
  ld <- as.matrix(fread(ld_file, header = FALSE))
  
  # Create SNP IDs in the format "chr:position" from the snp_names file
  snp_ids <- paste(chr, snp_names$V1, sep = ":")
  rownames(ld) <- snp_ids
  colnames(ld) <- snp_ids
  
  # Ensure the LD matrix is symmetric
  if (!isSymmetric(ld)) {
    cat("LD matrix is not symmetric; enforcing symmetry by averaging the lower and upper triangles.\n")
    ld <- (ld + t(ld)) / 2
  }
  
  # Identify common variants between the LD matrix and dt_window
  common_vars <- intersect(rownames(ld), dt_window$variant_id)
  if (length(common_vars) == 0) {
    stop("No common variants found between LD matrix and dt_window!")
  } else {
    cat("Number of common variants found:", length(common_vars), "\n")
  }
  
  # Subset the LD matrix to only include the common variants (and enforce dropping any extra rows/columns)
  ld <- ld[common_vars, common_vars, drop = FALSE]
  
  # Match the slope and standard error values using the order in common_vars.
  slope <- dt_window$slope[match(common_vars, dt_window$variant_id)]
  se    <- dt_window$slope_se[match(common_vars, dt_window$variant_id)]
  
  # Check for any NA values in slope or se
  na_idx <- which(is.na(slope) | is.na(se))
  if (length(na_idx) > 0) {
    cat("Dropping", length(na_idx), "variant(s) due to NA slope or se.\n")
    # Remove these indices from common_vars, slope, and se
    common_vars <- common_vars[-na_idx]
    slope <- slope[-na_idx]
    se <- se[-na_idx]
    ld <- ld[common_vars, common_vars, drop = FALSE]
    dt_window <- dt_window[variant_id %in% common_vars]
  }
  
  # Check that the diagonal of the LD matrix has no NAs (should be 1 for a correlation matrix)
  if (any(is.na(diag(ld)))) {
    cat("Found NA values in the diagonal of the LD matrix. Replacing NA values with 1.\n")
    diag(ld)[is.na(diag(ld))] <- 1
  }
  
  cat("Number of variants after dropping NA rows/columns:", length(common_vars), "\n")
  
  # (Optional) Verify that the LD matrix is still symmetric and the diagonal is correct.
  if (!isSymmetric(ld)) {
    stop("LD matrix is not symmetric after subsetting!")
  }
  if (!all(diag(ld) == 1)) {
    cat("Warning: Diagonal values of the LD matrix are not all 1. Inspecting...\n")
    print(diag(ld))
  }
  
  # Finally, run SuSiE with the LD matrix
  susie_fit <- susieR::susie_rss(L = 10, 
                                 iterations = 1000,
                                 estimate_residual_variance = T,
                                 z = slope/se,
                                 R = ld, 
                                 n = N1)
  
  # Create a copy of the LD matrix and set the lower triangle and diagonal to NA.
  #ld_upper <- ld
  #ld_upper[!upper.tri(ld_upper)] <- NA  # This sets the lower triangle AND the diagonal to NA
  #ld_long <- melt(ld_upper, varnames = c("SNP1", "SNP2"), value.name = "LD", na.rm = TRUE)
  # Plot the upper triangle of the LD matrix with aspect ratio 1
  # p1 <- {ggplot(ld_long, aes(x = SNP1, y = SNP2, fill = LD^2)) +
  #  geom_raster() +
  #  scale_fill_gradient2(midpoint = 0.5, low = "blue", mid = "white", high = "red", na.value = "grey50", name = "LD\nValue") +
  #  labs(title = "LD Matrix Heatmap", x = "SNP", y = "SNP") +
  #  coord_fixed() + theme_minimal() +
  #  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank())}
  
  # susie_plot(susie_fit, y = "PIP",  add_legend = T)
  # QC - check consistency of zscore and LD matrix
  zscore <- slope / se
  est <- estimate_s_rss(z = zscore, R = ld, n = N1)
  # condz_in = kriging_rss(z = zscore, R = ld, n = N1); condz_in$plot
  
  # Conditional get CS
  if (length(susie_fit$sets$cs) == 0 || all(sapply(susie_fit$sets$cs, length) == 0)) {
    rs <- data.table(gene = focalGene, 
                     SNP = NA, 
                     PIP = NA, 
                     estimate_s = est, 
                     nSNPs = as.numeric(length(common_vars)))
  } else {
    rs <- data.table(gene = focalGene,
                     SNP = rownames(ld_matrix)[unlist(susie_fit$sets$cs)],
                     PIP = susie_fit$pip[unlist(susie_fit$sets$cs)],
                     estimate_s = est,
                     nSNPs = as.numeric(length(common_vars)))
  }
  
  print(rs)
  return(rs)
}

# Run Susie
susie_results <- map(snps$common_gene, possibly(function(genei) {
  susieWindow(dt1 = qtl, focalGene = genei, N1 = N1)}, 
  otherwise = NULL))

susie_results_valid <- keep(susie_results, is.data.frame)
fin_susie <- rbindlist(susie_results_valid, fill = TRUE)

fwrite(fin_susie, file = out_file, sep = "\t")
