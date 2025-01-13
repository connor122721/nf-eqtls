#!/usr/bin/env Rscript
# Converts VCF files to a GDS file using SeqArray

# Libraries
library(SeqArray)
library(parallel)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check for correct number of arguments
if(length(args) != 3){
  stop("Usage: vcf_to_gds.R <input_vcf_files> <output_gds> <threads>")
}

input_vcf <- strsplit(args[1], ",")[[1]] # Comma-separated list
output_gds <- args[2]
threads <- as.numeric(args[3])

# Register cores
seqParallelSetup(threads)

# Convert VCF to GDS
seqVCF2GDS(vcf.fn = input_vcf, out.fn = output_gds, parallel = threads, verbose = TRUE)