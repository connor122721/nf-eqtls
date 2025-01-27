# Connor Murray
# Started 11.18.2024, edited 1.27.2025
# analyzing TOPchef eQTLs and perform colocalization
# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; R

# Libraries
library(data.table)
library(coloc)
library(tidyverse)
library(foreach)

# Argument parser
parser <- ArgumentParser()
parser$add_argument("--gwas", required=TRUE, help="Processed GWAS summary statistics.")
parser$add_argument("--eqtl", required=TRUE, help="cis-eQTL statistics.")
parser$add_argument("--chromosome", required=TRUE, help="Current chromosome.")
parser$add_argument("--N_gwas", required=TRUE, help="Samples in GWAS.")
parser$add_argument("--N_eqtl", required=TRUE, help="Samples in eQTL study.")

args <- parser$parse_args()

gwas_input <- args$gwas
eqtl_inpt <- args$eqtl


gwas_input="data/processed_GWAS/processed_levin22_gwas_HF_chr13.rds"
eqtl_inpt="nextflow_dna/output/tensorqtl_nominal/topchef_chr13_MaxPC30.cis_qtl_pairs.chr13.parquet"

### Datasets & Setup ###

# By chromosome coloc
coloc_chrom <- function(chr) {
  
  # Start on chromosome: chr="chr6"
  chromi=chr
  print(chromi)

  # Read in GWAS 
  gwas <- data.table(readRDS(gwas_input))
  
  # Restrict to all GWAS SNPs
  gwasSNPS <- unique(gwas$snpID_hg38)
  
  # Read in eQTL cis-nominal results (dense sampling)
  qtl <- data.table(arrow::read_parquet(eqtl_inpt))
  
  # Filter for significance and duplicated annotations
  qtl <- data.table(qtl[!duplicated(variant_id)][variant_id %in% gwasSNPS] %>% 
            mutate(chrom = chromi,
                   snp = as.numeric(tstrsplit(variant_id, ":")[[2]]),
                   type = "quant",
                   maf = case_when(af > 0.5 ~ 1-af,
                                   TRUE ~ af),
                   n = round(ma_count/(2*maf))))
  
  # Ensure same length
  gwas <- gwas[snpID_hg38 %in% qtl$variant_id]

  # Progress
  print(paste("Done reading in data:", chromi))
  
  ##### Make function to run coloc #####
  
  # Make list for colocalization (TOPChef)
  d1 <- list(beta=qtl$beta,
             varbeta=(qtl$beta_se^2),
             pvalues=qtl$pval,
             MAF=as.numeric(qtl$maf),
             N=N_eqtl,
             snp=qtl$variant_id,
             chrom=qtl$chrom,
             position=qtl$snp,
             type="quant")
  
  # Make list for colocalization (Levin)
  d2 <- list(beta=gwas$beta,
             varbeta=(gwas$standard_error^2),
             pvalues=as.numeric(gwas$p_value),
             MAF=as.numeric(gwas$maf),
             N=N_gwas,
             snp=gwas$snpID_hg38,
             chrom=gwas$chr_hg38,
             position=gwas$pos_hg38,
             type="quant")
  
  # Window size coloc function between Levin and eQTL datasets
  colocWindow <- function(dt1, N1=N_eqtl, dt2, N2=N_gwas, focalSNP, window=1e6) {
    # dt1=qtl; dt2=lewin; focalSNP="chr6:79027389[b37]"; window=1000000; N1=502; N2=1665481
  
    # Get focal SNP information
    pos = as.numeric(tstrsplit(focalSNP, ":")[[2]])
    Min_pos = pos-window
    Max_pos = pos+window
    
    # Filter datasets to focal SNPs
    dt1 <- dt1[chrom == chromi & snp >= Min_pos & snp <= Max_pos]
    dt2 <- dt2[chr_hg38 == chromi & pos_hg38 >= Min_pos & pos_hg38 <= Max_pos]
    focalGene = unique(dt1[variant_id==focalSNP]$phenotype_id)
    
    # Make list for colocalization (TOPChef)
    dt1 <- list(beta = dt1$slope,
                varbeta = dt1$slope_se^2,
                snp = dt1$variant_id,
                position = dt1$snp,
                type = "quant",
                N = N1,
                MAF = as.numeric(dt1$maf),
                pvalues = as.numeric(dt1$pval_nominal))
    
    # Make list for colocalization (Levin)
    dt2 <- list(beta = dt2$beta,
                varbeta = dt2$standard_error^2,
                snp = dt2$snpID_hg38,
                position = dt2$pos_hg38,
                type = "quant",
                N = N2,
                MAF = as.numeric(dt2$maf),
                pvalues = dt2$p_value)
    
    # Colocalization analysis using the coloc package
    coloc_res <- coloc.abf(dataset1 = dt1,
                           dataset2 = dt2)

    # Results of colocalization
    co <- data.table(chrom=chromi,
                     focalSNP=focalSNP,
                     gene=focalGene,
                     min_p.eqtl=min(dt1$pvalues),
                     min_p.gwas=min(dt2$pvalues),
                     nsnps=coloc_res$summary["nsnps"],
                     PP.H0=coloc_res$summary["PP.H0.abf"],
                     PP.H1=coloc_res$summary["PP.H1.abf"],
                     PP.H2=coloc_res$summary["PP.H2.abf"],
                     PP.H3=coloc_res$summary["PP.H3.abf"],
                     PP.H4=coloc_res$summary["PP.H4.abf"])
    
    # Finish
    return(co)
  }
  
  # Run coloc analysis with error handling
  coloc_results <- map(snps$V1, possibly(function(snp) {
    colocWindow(dt1 = qtl, 
                dt2 = gwas, 
                focalSNP = snp)}, otherwise = NULL))
  
  # Filter out NULL results and keep only data frames
  coloc_results_valid <- keep(coloc_results, is.data.frame)
  
  # Combine results with rbindlist
  fin_coloc <- rbindlist(coloc_results_valid, fill = TRUE)
  
  # Finish
  print(paste("Finish:", chr))
  return(fin_coloc)

}

# Run code by chromosome
dt = coloc_chrom(chromosome)

# Output results
write_delim(dt, file = paste("coloc_eqtl_levin22_HF_seVar", chromosome, ".txt", sep=""), delim = "\t")
