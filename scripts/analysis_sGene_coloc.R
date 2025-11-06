# Connor Murray
# Started 12.10.2024; modifed 1.30.2025
# analyzing TOPchef e/sQTLs and make coloc plots
# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

# Libraries
library(data.table)
library(tidyverse)
library(argparse)

# Argument parser
parser <- ArgumentParser()
parser$add_argument("--gwas", required=TRUE, help="Processed GWAS summary statistics, uses prefix to find relevant GWAS.")
parser$add_argument("--eqtl", required=TRUE, help="Number of max RNA PCs used in eQTL model.")
parser$add_argument("--chromosome", required=TRUE, help="Current chromosome.")
parser$add_argument("--gene", required=TRUE, help="Current gene.")

# Create list for LD matrix generation
# coloc <- fread("/scratch/csm6hg/nextflow_dna/output/coloc_sqtl/coloc_sqtl_candidates_full.txt") %>% 
#  filter(PP.H4 > 0.8) %>% select(chrom, minPos, maxPos, gene_raw, gene.x)
# write.table(coloc, file="/standard/vol185/cphg_Manichaikul/users/csm6hg/data/candgenes.spliceldlist_9_16_25.txt", sep = "\t", quote = F, row.names = F)

args <- parser$parse_args()

gwas_input <- as.character(args$gwas)
eqtl_inpt <- as.character(args$eqtl)
chromosome <- as.character(args$chromosome)
genei <- as.character(args$gene)

#setwd("/standard/vol185/cphg_Manichaikul/users/csm6hg/nextflow_dna/");gwas_input="shah20_gwas_HF";eqtl_inpt="MaxPC49";chromosome="chr10";genei="BAG3"
#setwd("/standard/vol185/cphg_Manichaikul/users/csm6hg/nextflow_dna/");gwas_input="jurgens24_gwas_HF";eqtl_inpt="MaxPC70";chromosome="chr7";genei="FLNC"
#setwd("/scratch/csm6hg/nextflow_dna/");gwas_input="jurgens24_gwas_HF";eqtl_inpt="MaxPC10";chromosome="chr2";genei="chr2:178678497:178678747:clu_58053_-:ENSG00000155657.27"

# Make sure the GWAS input is correct
if (gwas_input == "jurgens24_gwas_HF") {
    gwasName <- "Jurgens 2024"
}

print(gwas_input)
print(eqtl_inpt)
print(chromosome)
print(genei)

### Datasets & Setup ###
gtf <- data.table(readRDS("/standard/vol185/cphg_Manichaikul/users/csm6hg/genome_files/gencode.v34.GRCh38.ERCC.genes.collapsed.streamlined.RDS")) %>% 
  select(-c(V7:V9, file), gene_length=len)

# Common theme element
themei <- {
  theme_bw() + 
    theme(axis.title.x = element_text(face = "bold", size = 18),
          axis.text.x = element_text(face = "bold", size = 16),
          axis.title.y = element_text(face = "bold", size = 18),
          axis.text.y = element_text(face = "bold", size = 16),
          legend.text = element_text(face = "bold", size = 16),
          legend.title = element_text(face = "bold", size = 18)) 
}

#### MAKE LD aware plots ####
print("reading in data")

# Read in topchef QTLs and restrict to colocalized genes
qtl <- data.table(arrow::read_parquet(list.files(path = "/scratch/csm6hg/nextflow_dna/output/splicing/sqtls_nominal/", 
                  pattern = paste(chromosome, eqtl_inpt, sep="_"), full.names = T)))
prefix="sqtl"

# Read in GWAS Shah 2020
gwas <- data.table(readRDS(list.files(path = "/scratch/csm6hg/nextflow_dna/output/gwas/", 
                   pattern = paste(gwas_input, "_", chromosome, ".rds", sep=""), full.names = T)))

# Merge
dt <- data.table(na.omit(qtl %>%
        filter(phenotype_id==genei) %>% 
        left_join(gwas, 
              by=c("variant_id"="snpID_hg38")) %>% 
        mutate(gene=tstrsplit(phenotype_id, ":")[[5]],
               gene_edit=tstrsplit(gene,".", fixed=T)[[1]]) %>% 
        left_join(gtf, by=c("gene_edit", "gene"))))

print("finished reading in data")

# Read in LD matrix for candidate gene!
ld.file <- list.files(path = "/scratch/csm6hg/nextflow_dna/output/splice_linkage", 
                      pattern = genei, 
                      full.names = T)

# Read and process LD file
readLD_process <- function(ldFile) {
  
  # Read in LD file: ldFile=ld.file
  if (length(ldFile) > 1) {
    ldFile <- ldFile[grepl("ld\\.gz$", ldFile)]
    if (length(ldFile) == 0) {
      stop("No gzipped LD file (ld.gz) found among the provided files!")
    }
    ldFile <- ldFile[1]  # if multiple, take the first one
  }
  
  snp_names <- fread(ld.file[ld.file %like% "snp_names.txt"], header = FALSE)
  
  # Read in the LD matrix (fread can handle gzipped files)
  ld <- as.matrix(fread(ldFile, header = FALSE))
  
  # Create SNP IDs in the format "chr:position" from the snp_names file
  snp_ids <- paste(chromosome, snp_names$V1, sep = ":")
  rownames(ld) <- snp_ids
  colnames(ld) <- snp_ids
  
  # Ensure the LD matrix is symmetric
  if (!isSymmetric(ld)) {
    cat("LD matrix is not symmetric; enforcing symmetry by averaging the lower and upper triangles.\n")
    ld <- (ld + t(ld)) / 2
  }
  
  # Identify common variants between the LD matrix and dt_window
  common_vars <- intersect(rownames(ld), dt[phenotype_id==genei]$variant_id)
  if (length(common_vars) == 0) {
    stop("No common variants found between LD matrix and dt_window!")
  } else {
    cat("Number of common variants found:", length(common_vars), "\n")
  }
  
  # Subset the LD matrix to only include the common variants (and enforce dropping any extra rows/columns)
  ld <- ld[common_vars, common_vars, drop = FALSE]
  
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
  
  # Convert the LD matrix to a data.table
  ld_dt <- as.data.table(as.table(ld))
  setnames(ld_dt, c("V1", "V2", "N"), c("SNPA", "SNPB", "R"))
  ld_dt <- ld_dt %>% mutate(R2=R^2)

  # Finish
  return(ld_dt)
}

# Process
ld_dt <- data.table(readLD_process(ld.file))

# Find candidate sentinenal variant
dt_sub <- dt[phenotype_id %in% genei]
dt_sub[, rank_eqtl := rank(pval_nominal), by = common_gene]
dt_sub[, rank_gwas := rank(p_value), by = common_gene]
dt_sub[, total_rank := rank_eqtl + rank_gwas]

# For each gene, choose the variant with the smallest total_rank
sentinel_combined <- dt_sub[order(rank_gwas), .SD[1], by = common_gene]
sent = sentinel_combined$variant_id

# Helper function to make HTTP GET requests
library(httr)
make_request <- function(url) {
  response <- GET(url)
  if (status_code(response) != 200) {
    stop(sprintf("Request to %s failed with status code: %s", url, status_code(response)))
  }
  response
}

# Function to retrieve the rsID(s) given a chromosome and position
get_rsid <- function(chrom, pos) {
  
  # Construct the variation URL (GFF3 format)
  variation_url <- sprintf("https://rest.ensembl.org/overlap/region/homo_sapiens/%s:%s-%s?feature=variation;content-type=text/x-gff3",
                           chrom, pos, pos)
  
  # Make the API request
  response <- make_request(variation_url)
  variation_text <- content(response, as = "text", encoding = "UTF-8")
  
  # Initialize vector to store rsIDs
  rsid_list <- c()
  
  # Process the response line by line
  lines <- unlist(strsplit(variation_text, "\n"))
  for (line in lines) {
    if (nchar(line) == 0 || startsWith(line, "#")) next  # Skip empty or comment lines
    
    # Split the line by tabs and check if there are enough fields
    fields <- strsplit(line, "\t")[[1]]
    if (length(fields) >= 9) {
      attributes <- fields[9]
      
      # Split the attribute field by ";" and search for the ID attribute
      attrs <- unlist(strsplit(attributes, ";"))
      for (attr in attrs) {
        attr <- trimws(attr)
        if (startsWith(attr, "ID=")) {
          rsid <- sub("ID=", "", attr)
          # Remove any "sequence_variant:" prefix if present
          rsid <- sub("sequence_variant:", "", rsid)
          rsid_list <- c(rsid_list, rsid)
          break  # Stop after the first ID attribute in this line
        }
      }
    }
  }
  
  # Return unique rsIDs (in case duplicates occur)
  unique(rsid_list)
}

# Example usage:
rsid <- get_rsid(sentinel_combined$chrom, sentinel_combined$pos_hg38)
print(rsid)

# Merge with gene
dt1 <- dt[phenotype_id %in% genei] %>% 
  left_join(ld_dt %>% 
              filter(SNPA == sent), 
            by=c("variant_id"="SNPB"))

if( length(rsid) >= 1  ) {
  dt1$lead_rsid <- rsid[1]
}

if( length(rsid) == 0 ){
  dt1$lead_rsid <- sent
}

# Plot gwas and eqtl statistics
library(ggrepel)
library(patchwork)

# Make LocusComparer plot
plot1 <- {
  
  dt1 %>% 
    ggplot(aes(x = -log10(pval_nominal),
               y = -log(p_value),
               color = R2)) +
    geom_point(size=2, alpha=0.6) +
    geom_point(data = dt1[variant_id == sent], 
               size = 5, shape = 8, color = "black") +
    geom_label_repel(data = dt1[variant_id == sent], 
                     aes(label = lead_rsid), alpha = 0.8) +
    facet_wrap(~common_gene, scales = "free") +
    scale_color_gradient2(low = "blue", mid = "green", high = "red", 
                          midpoint = 0.5, name = expression(italic(R)^2)) +
    labs(x = expression("-log(cis-eQTL " * italic(p) * "-value)"),
         y = expression("-log(GWAS " * italic(p) * "-value)")) +
    themei +
    theme(aspect.ratio = 1,
          strip.text = element_text(face = "bold.italic", size = 20),
          legend.title = element_text(face = "bold", size = 20))

}

# eQTL Manhatten plot
plot_eqtl <- {
  
  dt1 %>% 
    ggplot(aes(x = pos_hg38/1e6,
               y = -log(pval_nominal),
               color = R2)) +
    geom_point(size=2, alpha=0.6) +
    geom_point(data = dt1[variant_id == sent], 
               size = 5, shape = 8, color = "black") +
    geom_label_repel(data = dt1[variant_id == sent], 
                     aes(label = lead_rsid), alpha = 0.8) +
    facet_wrap(~common_gene, scales = "free") +
    scale_color_gradient2(low = "blue", mid = "green", high = "red", 
                          midpoint = 0.5, name = expression(italic(R)^2)) +
    labs(x = "Position (Mbps)",
         y = expression("-log(TOPCHef eQTL " * italic(p) * "-value)")) +
    themei +
    theme(aspect.ratio = 1,
          strip.text = element_text(face = "bold.italic", size = 20),
          legend.title = element_text(face = "bold", size = 20))
  
}

# GWAS Manhatten plot
plot_gwas <- {
  
  dt1 %>% 
    ggplot(aes(x = pos_hg38/1e6,
               y = -log(p_value),
               color = R2)) +
    geom_point(size=2, alpha=0.6) +
    geom_point(data = dt1[variant_id == sent], 
               size = 5, shape = 8, color = "black") +
    geom_label_repel(data = dt1[variant_id == sent], 
                     aes(label = lead_rsid), alpha = 0.8) +
    facet_wrap(~common_gene, scales = "free") +
    scale_color_gradient2(low = "blue", mid = "green", high = "red", 
                          midpoint = 0.5, name = expression(italic(R)^2)) +
    labs(x = "Position (Mbps)",
         y = expression("-log(GWAS " * italic(p) * "-value)")) +
    themei +
    theme(aspect.ratio = 1,
          strip.text = element_text(face = "bold.italic", size = 20),
          legend.title = element_text(face = "bold", size = 20))
  
}

# Gene regions subplot!
library(ggrepel)

gtf.dt <- gtf[chrom==unique(dt1$chromosome)][start %in% min(dt1$pos_hg38):max(dt1$pos_hg38)][common_gene %in% unique(dt1$common_gene)]

gene_reg <- {
  
  gtf.dt %>% 
    ggplot(aes(x = start/1e6, xend = stop/1e6,
               y = reorder(common_gene, start),
               yend = reorder(common_gene, start))) +
    geom_segment(color = "lightblue2", linewidth = 4) +
    geom_text_repel(data = gtf.dt[start %in% (unique(dt1$start)-200000):(unique(dt1$start)+200000)], 
                    aes(x = (start + stop) / (2 * 1e6), label = common_gene),
              color = "black", size = 3, vjust = -0.5) +
    geom_vline(xintercept = dt1[SNP==lead_rsid]$pos_hg38/1e6, linetype=4, color="red") +
    themei +    
    theme_classic() +
    theme(aspect.ratio = 0.3, 
          axis.title.x = element_text(face = "bold", size = 18),
          axis.text.x = element_text(face = "bold", size = 16),
          axis.title.y = element_text(face = "bold", size = 18),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    labs(x = "Position (Mbps)", y = "Gene")

}

combined_plot <- (plot1 / plot_eqtl / plot_gwas) / gene_reg

# Save Output!
ggsave(plot = combined_plot, 
       filename = paste(genei, ".topchef.coloc.", gwas_input, ".pdf", sep=""), 
       width = 20, height = 20)
