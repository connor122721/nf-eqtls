#!/usr/bin/env Rscript
# Converts VCF files to a GDS file using SeqArray

# Libraries
library(SeqArray)
library(parallel)
library(argparse)

# Argument parser
parser <- ArgumentParser()
parser$add_argument("--input_vcf_files", required=TRUE, help="Comma-separated list of input VCF files")
parser$add_argument("--output_gds", required=TRUE, help="Path to output GDS file")
parser$add_argument("--threads", required=TRUE, type="integer", help="Number of threads to use")
args <- parser$parse_args()

input_vcf <- strsplit(args$input_vcf_files, ",")[[1]]
output_gds <- args$output_gds
threads <- args$threads

# Register cores
seqParallelSetup(threads)

# Convert VCF to GDS
seqVCF2GDS(vcf.fn = input_vcf, out.fn = output_gds, parallel = threads, verbose = TRUE)