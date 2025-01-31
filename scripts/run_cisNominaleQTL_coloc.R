# Connor Murray
# Started 12.8.2024; modified 1.30.2025
# analyzing TOPchef eQTLs and perform colocalization w/HF GWAS
# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

# Libraries
library(data.table)
library(coloc)
library(tidyverse)
library(foreach)
library(argparse)

# Argument parser
parser <- ArgumentParser()
parser$add_argument("--gwas", required=TRUE, help="Processed GWAS summary statistics.")
parser$add_argument("--eqtl", required=TRUE, help="cis-eQTL statistics.")
parser$add_argument("--shortList", required=TRUE, help="Extracts candidate SNPs/eGenes to test for colocalization from first cis eqtl-tensorQTL run.")
parser$add_argument("--chromosome", required=TRUE, help="Current chromosome.")
parser$add_argument("--N_gwas", required=TRUE, help="Samples in GWAS.")
parser$add_argument("--N_eqtl", required=TRUE, help="Samples in eQTL study.")
parser$add_argument("--prefix", required=TRUE, help="Prefix of output files.")

args <- parser$parse_args()

gwas_input <- args$gwas
eqtl_inpt <- args$eqtl
sig_genes <- unlist(fread(args$shortList) %>% 
                      select(phenotype_id))
chromosome <- args$chromosome
N_gwas <- as.numeric(args$N_gwas)
N_eqtl <- as.numeric(args$N_eqtl)
output_pre <- args$prefix

#gwas_input="output/gwas/processed_levin22_gwas_HF_chr10.rds"
#eqtl_inpt="output/tensorqtl_nominal/topchef_chr10_MaxPC49.cis_qtl_pairs.chr10.parquet"
#shortlist="output/tensorqtl/topchef_chr10_MaxPC49.cis_qtl.txt.gz"
#N_gwas=1665481
#N_eqtl=516
#chromosome="chr10"

### Datasets & Setup ###

# By chromosome coloc
coloc_chrom <- function(chr) {
  
  # Start on chromosome: chr="chr10"
  chromi=chr
  print(paste("Running:", chromi, sep=" "))
  
  # Read in parquet files
  qtl <- data.table(arrow::read_parquet(eqtl_inpt)) %>% 
    mutate(chrom = chromi,
           snp = as.numeric(str_remove_all(tstrsplit(variant_id, ":")[[2]], "\\[b37\\]")),
           type = "quant",
           maf = case_when(af > 0.5 ~ 1-af,
                           TRUE ~ af),
           n = round(ma_count/(2*maf)))
  
  # Read in Levin et al 2022 GWAS
  gwas <- data.table(readRDS(gwas_input))
  
  # Progress
  print(paste("Done reading/preparing data:", chromi))
  
    # Window size coloc function between Shah and eQTL datasets
    colocWindow <- function(dt1, N1=N_eqtl, focalGene, dt2, N2=N_gwas) {
      # dt1=qtl; dt2=gwas; focalGene="ENSG00000177791"; N1=516; N2=1665481
      
      # Restrict to all eQTL SNPs 
      dt1 <- dt1[phenotype_id %in% focalGene][!duplicated(variant_id)]
      min_window <- as.numeric(min(dt1$snp))
      max_window <- as.numeric(max(dt1$snp))
      
      # Ensure same length
      dt2 <- dt2[snpID_hg38 %in% dt1$variant_id]
      dt1 <- dt1[variant_id %in% dt2$snpID_hg38]
      
      # Make list for colocalization (TOPChef)
      dt1 <- list(snp = dt1$variant_id,
                  position = dt1$snp,
                  type = "quant",
                  N = N1,
                  MAF = as.numeric(dt1$maf),
                  pvalues = as.numeric(dt1$pval_nominal))
      
      # Make list for colocalization (Shah)
      dt2 <- list(snp = dt2$snpID_hg38,
                  position = dt2$pos_hg38,
                  type = "quant",
                  N = N2,
                  MAF = as.numeric(dt2$maf),
                  pvalues = as.numeric(dt2$p_value))
      
      # Colocalization analysis using the coloc package
      coloc_res <- coloc.abf(dataset1 = dt1,
                             dataset2 = dt2)
      
      maxi = max(coloc_res$results$SNP.PP.H4)
      snpi = coloc_res$results[which(coloc_res$results$SNP.PP.H4==maxi),]$snp
      
      #plot(coloc_res)

      # Results of colocalization
      co <- data.table(chrom=chromi,
                       minPos=min_window,
                       maxPos=max_window,
                       gene=focalGene,
                       min_p.eqtl=min(dt1$pvalues),
                       min_p.gwas=min(dt2$pvalues),
                       nsnps=coloc_res$summary["nsnps"],
                       PP.H0=coloc_res$summary["PP.H0.abf"],
                       PP.H1=coloc_res$summary["PP.H1.abf"],
                       PP.H2=coloc_res$summary["PP.H2.abf"],
                       PP.H3=coloc_res$summary["PP.H3.abf"],
                       PP.H4=coloc_res$summary["PP.H4.abf"],
                       H4_H3_ratio=coloc_res$summary["PP.H4.abf"]/coloc_res$summary["PP.H3.abf"],
                       maxSNP=snpi,
                       maxPP.H4=maxi,
                       gwas=output_pre)
      
      # Finish
      return(co)
    }
    
  # Run coloc analysis with error handling
  coloc_results <- map(sig_genes, possibly(function(gene) {
    colocWindow(dt1 = qtl, 
                dt2 = gwas, 
                focalGene = gene)}, otherwise = NULL))
  
  # Filter out NULL results and keep only data frames
  coloc_results_valid <- keep(coloc_results, is.data.frame)
  
  # Combine results with rbindlist
  fin_coloc <- rbindlist(coloc_results_valid, fill = TRUE)
  
  # Finish
  print(paste("Finish:", chr))
  return(fin_coloc)

}

# Run code
dt <- coloc_chrom(chromosome)

# Output results
write_delim(dt, file = paste("coloc_eqtl_", output_pre, "_", chromosome, ".txt", sep=""), delim = "\t")
