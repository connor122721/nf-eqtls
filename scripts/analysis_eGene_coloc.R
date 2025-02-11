# Connor Murray
# Started 12.10.2024; modifed 1.30.2025
# analyzing TOPchef eQTLs and make 
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

args <- parser$parse_args()

gwas_input <- as.character(args$gwas)
eqtl_inpt <- as.character(args$eqtl)
chromosome <- as.character(args$chromosome)
genei <- as.character(args$gene)

#setwd("/standard/vol185/cphg_Manichaikul/users/csm6hg/nextflow_dna/");gwas_input="shah20_gwas_HF";eqtl_inpt="MaxPC49";chromosome="chr10";genei="BAG3"

# Make sure the GWAS input is correct
if (gwas_input == "jurgens24_gwas_HF") {
    gwasName <- "Jurgens 2024"
} else if (gwas_input == "levin22_gwas_HF") {
    gwasName <- "Levin 2022"
} else if (gwas_input == "shah20_gwas_HF") {
    gwasName <- "Shah 2020"
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

# Read in topchef eQTLs and restrict to colocalized genes
qtl <- data.table(arrow::read_parquet(list.files(path = "/standard/vol185/cphg_Manichaikul/users/csm6hg/nextflow_dna/output/tensorqtl_nominal/", 
                                                 pattern = paste(chromosome, eqtl_inpt, sep="_"), full.names = T)))

# Read in GWAS Shah 2020
gwas <- data.table(readRDS(list.files(path = "/standard/vol185/cphg_Manichaikul/users/csm6hg/nextflow_dna/output/gwas/", 
                                      pattern = paste(gwas_input, "_", chromosome, ".rds", sep=""), full.names = T)))

# Merge
dt <- data.table(na.omit(qtl %>%
        left_join(gwas, 
              by=c("variant_id"="snpID_hg38"))) %>% 
        left_join(gtf, by=c("phenotype_id"="gene_edit")))

print("finished reading in data")

# Read in LD matrix for candidate gene!
ld.file <- list.files(path = "/standard/vol185/cphg_Manichaikul/users/csm6hg/nextflow_dna/output/linkage", 
                      pattern = genei, 
                      full.names = T)

# Read and process LD file
readLD_process <- function(ldFile) {
  
  # Read in LD file: ldFile=ld.file
  ld <- data.table(fread(input = ldFile, header = T) %>% 
                     mutate(SNPA = paste("chr", CHR_A, ":", BP_A, sep=""),
                            SNPB = paste("chr", CHR_B, ":", BP_B, sep="")) %>% 
                     filter(SNPA %in% dt$variant_id & SNPB %in% dt$variant_id) %>% 
                     mutate(R = sqrt(R2)) %>% 
                     select(SNPA, SNPB, R, R2))
  
  # Get all unique SNPs
  snp_list <- unique(c(ld$SNPA, ld$SNPB))
  snp_list <- sort(snp_list) # Sort for consistency
  
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

# Find candidate sentinenal variant
dt_sub <- dt[common_gene %in% genei]
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
dt1 <- na.omit(dt[common_gene %in% genei] %>% 
  left_join(ld_dt %>% 
              filter(SNPA == sent), 
            by=c("variant_id"="SNPB")))

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
    theme(strip.text = element_text(face = "bold.italic", size = 20),
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
    theme(strip.text = element_text(face = "bold.italic", size = 20),
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
    theme(strip.text = element_text(face = "bold.italic", size = 20),
          legend.title = element_text(face = "bold", size = 20))
  
}

combined_plot <- (plot1 / plot_eqtl / plot_gwas)

# Save Output!
ggsave(plot = combined_plot, 
       filename = paste(genei, ".topchef.coloc.", gwas_input, ".pdf", sep=""), 
       width = 7, height = 14)
