#!/usr/bin/env python
# coding: utf-8
"""
Author: Connor Murray
Date: 2024-12-14

Description:
  This script reads in a GTF file, extracts TSS positions for each gene, and
  writes them as a BED file for downstream eQTL workflows.
  
Usage:
  ./prep_gtf_for_tss.py --gtf <path_to_input_GTF> --output <path_to_output_BED>
"""
import argparse
import pandas as pd
import qtl.io as qtl_io
import qtl.torus as qtl_torus

def main():
    parser = argparse.ArgumentParser(
        description="Extract TSS from a GTF file and write to BED format.")
    parser.add_argument("--gtf", required=True, help="Path to GTF file.")
    parser.add_argument("--output", required=True, help="Path for output BED file.")
    args = parser.parse_args()

    gtf_file = args.gtf
    output_file = args.output

    # Print input info
    print(f"Reading GTF: {gtf_file}")

    # Optional: Show a small preview (comment out if large GTF)
    gtf_b = pd.read_csv(gtf_file, sep="\t", comment="#", header=None)
    print("Preview of raw GTF (first 5 lines):")
    print(gtf_b.head())

    # Convert GTF => TSS BED
    print("Converting GTF to TSS bed...")
    gtf_tss = qtl_io.gtf_to_tss_bed(
        gtf_file, 
        feature='gene', 
        phenotype_id='gene_id')

    print("Preview of TSS BED (first 5 lines):")
    print(gtf_tss.head())

    # Write out BED
    gtf_tss.to_csv(output_file, sep="\t", index=False)
    print(f"Saved TSS BED to: {output_file}")

if __name__ == "__main__":
    main()
