#!/usr/bin/env python
# coding: utf-8

"""
By: Connor S. Murray
  - Performs differential gene expression on raw RNA-seq data with DESeq2!
"""
# Load Libraries
import argparse
import sys
import os
import math
import pandas as pd
import numpy as np
import re
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from pydeseq2.utils import load_example_data
import matplotlib.pyplot as plt
import seaborn as sns

# Arguments from command line
def parse_arguments():
    parser = argparse.ArgumentParser(description="Differential gene expression analysis with DEseq2.")
    parser.add_argument("--metadata", required=False, 
                        default="/standard/vol185/cphg_Manichaikul/users/csm6hg/metadata/metadata_10_17_2024_CSM.txt",
                        help="Path to metadata file (must include columns: SAMPLE_ID_TOR, diagnosis_simple, Age_at_collection, Gender).")
    parser.add_argument("--gene_counts", required=False, 
                        default="/standard/vol185/TOPMed/TOPCHef/82214/topmed-dcc/exchange/phs002038_TOPMed_TOPCHeF/Omics/RNASeq/release3/TOPMed_Taylor_P4.RNASeQCv2.3.6_gene_reads.gct.gz",
                        help="Raw gene counts file")
    parser.add_argument("--gtf", required=False, 
                        default="/standard/vol185/cphg_Manichaikul/users/csm6hg/genome_files/gencode.v34.GRCh38.ERCC.genes.collapsed.gtf", 
                        help="Path to GTF bed file.")
    parser.add_argument("--outdir", required=False, 
                        default="/standard/vol185/cphg_Manichaikul/users/csm6hg/nextflow_dna/output/dge", 
                        help="Path to output directory")
    return parser.parse_args()

# TESTING 
#metadata="/standard/vol185/cphg_Manichaikul/users/csm6hg/metadata/metadata_10_17_2024_CSM.txt";gene_counts="/standard/vol185/TOPMed/TOPCHef/82214/topmed-dcc/exchange/phs002038_TOPMed_TOPCHeF/Omics/RNASeq/release3/TOPMed_Taylor_P4.RNASeQCv2.3.6_gene_reads.gct.gz";gtf="/standard/vol185/cphg_Manichaikul/users/csm6hg/genome_files/gencode.v34.GRCh38.ERCC.genes.collapsed.gtf"

# Function to extract gene_id from the attributes field using regex
def extract_gene_info(attr_str):
    m1 = re.search(r'gene_id "([^"]+)"', attr_str)
    m2 = re.search(r'gene_name "([^"]+)"', attr_str)
    gene_id = m1.group(1) if m1 else None
    gene_name = m2.group(1) if m2 else None
    # Remove trailing version information: a dot followed by one or more digits at the end
    gene = re.sub(r"\.[0-9]+$", "", gene_id) if gene_id else None
    return pd.Series([gene_id, gene, gene_name], index=["gene_id", "gene", "common_gene"])

# Define the column names for a standard GTF file
gtf_columns = ["chrom", "source", "feature", "start", "end", 
               "score", "strand", "frame", "attribute"]

# Read the GTF file
df = pd.read_csv(gtf, sep="\t", comment="#", header=None, names=gtf_columns)

# Optionally, filter the DataFrame to only include gene-level annotations
genes_df = df[df["feature"] == "gene"].copy()

# Extract gene_id from the attribute column
genes_df[["gene_id", "gene","common_gene"]] = genes_df["attribute"].apply(extract_gene_info)

# Select only the columns of interest
result = genes_df[["chrom", "start", "end", "gene_id", "gene", "common_gene"]]
print(result.head())

# Read in metadata
meta = pd.read_csv(metadata, sep="\t", index_col="SAMPLE_ID_TOR")
print(meta.head())

if "TOR238072" in meta.index:
        meta.drop("TOR238072", inplace=True)

# Read in raw gene count matrix
gene_count = pd.read_table(
    gene_counts,
    compression="gzip",
    skiprows=2,
    index_col="Name")
print(f"Original gene_count shape: {gene_count.shape}")

# Ensure columns are string
gene_count.columns = gene_count.columns.astype(str)

# Align with meta
common_samples = meta.index.intersection(gene_count.columns)
print(f"Number of common samples: {len(common_samples)}")

gene_count = gene_count.loc[:, common_samples]
print(f"Filt. gene_count shape: {gene_count.shape}")
print(gene_count.head())

# Extract pre-filtered gene list
genes_list_filt = pd.read_csv("../output/rna/norm_medrat.tsv", sep="\t")
genes_list = genes_list_filt.filter(like='ENSG').columns
gene_count = gene_count[gene_count.index.isin(genes_list)]
print(f"Filt. gene_count shape: {gene_count.shape}")

# Extract pre-filtered sample list
samples = pd.read_csv("../output/eqtl/topchef_samples_1_15_25.txt", header=None, names=["SAMPLE_ID_NWD"])

gene_countsi = gene_count.T
meta_re = meta.rename_axis('Name')

# Convert the SAMPLE_ID_NWD column to a list
sample_ids = samples["SAMPLE_ID_NWD"].tolist()

# Filter metadata rows where SAMPLE_ID_NWD is in the sample_ids list
meta_re = meta_re[meta_re["SAMPLE_ID_NWD"].isin(sample_ids)]
meta_re = meta_re[meta_re["diagnosis_simple"].isin(["IDCM", "Non-Failing"])] # Restrict to just DCM vs Control
design_columns = ['diagnosis_simple', 'Age_at_collection']
meta_re = meta_re.dropna(subset=design_columns)
gene_countsii = gene_countsi[gene_countsi.index.isin(meta_re.index)]

print(f"Gene count matrix for DESeq2: {gene_countsii.shape}")

# Make DDS object
inference = DefaultInference(n_cpus=4)
dds = DeseqDataSet(counts=gene_countsii, 
                   metadata=meta_re, 
                   design_factors=['diagnosis_simple', 'Age_at_collection', 'Gender'],
                   refit_cooks=True,
                   inference=inference,)

# Normalization
dds.fit_size_factors()
norm_counts = pd.DataFrame(dds.layers["normed_counts"])
norm_counts.index = gene_countsii.index
norm_counts.columns = gene_countsii.columns
print(f"Normalized count matrix shape: {norm_counts.shape}")

# Run differential expression analysis (this takes a couple of mins)
dds.deseq2()

# print(dds)
# print(meta["diagnosis_simple"].value_counts())

# Contrast groups 
stat_dcm = DeseqStats(dds, contrast = ('diagnosis-simple', 'IDCM', 'Non-Failing'), alpha=0.05)
#stat_icm = DeseqStats(dds, contrast = ('diagnosis-simple', 'ICM', 'Non-Failing'), alpha=0.05)
#stat_dcm.summary()

# PCA
import scanpy as sc
sc.tl.pca(dds)
sc.pl.pca(dds, 
          groups=["Non-Failing", "IDCM"], 
          color = 'diagnosis-simple', 
          ncols=2,
          size = 200, 
          components=['1,2', '2,3', '3,4', '4,5'])

# Extract significant genes
res = stat_dcm.results_df
res = res.merge(result, left_index=True, right_on="gene_id", how='left')
res = res[res.baseMean >= 10]
res.to_csv(os.path.join(outdir, 'sig_dge_DCM.csv'), index=True)
sigs = res[(res.padj < 0.05) & (abs(res.log2FoldChange) > math.log2(1.5))]
print(sigs.head())
