#!/usr/bin/env python3

# Does coloc on given candidate gene list and processed GWAS
# By: Connor S. Murray; 2.26.2025

# Import packages
import argparse
import pandas as pd
import numpy as np
import re
import pyreadr  # for reading RDS files
import pyarrow.parquet as pq  # alternatively, you can use pd.read_parquet
from qtl.coloc import abf  # assuming a package “qtl” with coloc_abf exists
from tqdm import tqdm
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

# Arguments from NF pipeline
def parse_args():
    parser = argparse.ArgumentParser(description="Analyze TOPchef eQTLs and perform colocalization with HF GWAS")
    parser.add_argument("--gwas", required=True, help="Processed GWAS summary statistics (RDS file).")
    parser.add_argument("--eqtl", required=True, help="cis-eQTL statistics (parquet file).")
    parser.add_argument("--shortList", required=True, help="Candidate SNPs/eGenes file from tensorQTL run (tab-delimited, possibly gzipped).")
    parser.add_argument("--chromosome", required=True, help="Current chromosome (e.g., 'chr7').")
    parser.add_argument("--N_gwas", required=True, type=int, help="Number of samples in GWAS.")
    parser.add_argument("--N_eqtl", required=True, type=int, help="Number of samples in eQTL study.")
    parser.add_argument("--prefix", required=True, help="Prefix of output files.")
    return parser.parse_args()

# TEST!
# gwas="output/gwas/processed_jurgens24_gwas_HF_chr10.rds";eqtl="output/tensorqtl_nominal/topchef_chr10_MaxPC49.cis_qtl_pairs.chr10.parquet";shortList="output/tensorqtl/topchef_chr10_MaxPC49.cis_qtl.txt.gz";chromosome="chr10"; N_gwas=955733; N_eqtl=516; prefix="test"

def extract_snp(variant_id):
    """Extract the SNP numeric position from variant_id string.
       Expected format: 'chr7:123456[b37]'."""
    try:
        parts = variant_id.split(":")
        if len(parts) < 2:
            return np.nan
        # Remove any bracketed annotations like "[b37]"
        number_str = re.sub(r'\[.*?\]', '', parts[1])
        return float(number_str)
    except Exception:
        return np.nan

def colocWindow_for_gene(focalGene, dt1, dt2, N1, N2, chrom, output_pre):
    """
    Process a single candidate gene.
    dt1: cis-eQTL DataFrame.
    dt2: GWAS DataFrame.
    """
    # Filter to the gene of interest and remove duplicate SNPs
    dt1_gene = dt1[dt1['phenotype_id'] == focalGene].drop_duplicates(subset=['variant_id'])
    print(f"Processing gene: {focalGene}")
    if dt1_gene.empty:
        return None
    min_window = dt1_gene['snp'].min()
    max_window = dt1_gene['snp'].max()
    
    # Restrict dt2 to SNPs found in dt1_gene and vice versa
    dt2_subset = dt2[dt2['snpID_hg38'].isin(dt1_gene['variant_id'])]
    dt1_gene = dt1_gene[dt1_gene['variant_id'].isin(dt2_subset['snpID_hg38'])]
    if dt1_gene.empty or dt2_subset.empty:
        return None
    
    # Prepare dataset for coloc analysis (eQTLs)
    dataset1 = {
        "snp": dt1_gene["variant_id"].tolist(),
        "position": dt1_gene["snp"].tolist(),
        "type": "quant",
        "sample_size": int(N1),
        "maf": dt1_gene["maf"].tolist(),
        "pval_nominal": dt1_gene["pval_nominal"].tolist(),
        "beta": dt1_gene["slope"].tolist(),
        "beta_se": dt1_gene["slope_se"].tolist()
    }
    
    # Prepare dataset for coloc analysis (GWAS data)
    dataset2 = {
        "snp": dt2_subset["snpID_hg38"].tolist(),
        "position": dt2_subset["pos_hg38"].tolist(),
        "type": "quant",
        "sample_size": int(N2),
        "maf": dt2_subset["maf"].tolist(),
        "pval_nominal": dt2_subset["p_value"].tolist(),
        "beta": dt2_subset["beta"].tolist(),
        "beta_se": dt2_subset["standard_error"].tolist()
    }
    
    d1 = pd.DataFrame(dataset1)
    d2 = pd.DataFrame(dataset2)
    
    # Run colocalization analysis and unpack the results
    summary, results_df = abf(df1=d1, df2=d2)
    
    # Extract the SNP with the maximum PP.H4
    maxi = results_df["snp_pp_h4"].max()
    snpi = results_df.loc[results_df["snp_pp_h4"] == maxi, "snp_1"].iloc[0]
    
    co = {
        "chrom": chrom,
        "minPos": min_window,
        "maxPos": max_window,
        "gene": focalGene,
        "min_p.eqtl": min(dataset1["pval_nominal"]),
        "min_p.gwas": min(dataset2["pval_nominal"]),
        "nsnps": d1.shape[0],
        "PP.H0": summary["pp_h0_abf"],
        "PP.H1": summary["pp_h1_abf"],
        "PP.H2": summary["pp_h2_abf"],
        "PP.H3": summary["pp_h3_abf"],
        "PP.H4": summary["pp_h4_abf"],
        "H4_H3_ratio": summary["pp_h4_abf"] / summary["pp_h3_abf"] if summary["pp_h3_abf"] != 0 else np.nan,
        "maxSNP": snpi,
        "maxPP.H4": maxi,
        "gwas_pre": output_pre
    }
    return co

def coloc_chrom(eqtl_inpt, gwas_input, sig_genes, chrom, N_eqtl, N_gwas, output_pre):
    
    print(f"Running: {chrom}")
    
    # Read in the cis-eQTL data from a parquet file
    qtl = pd.read_parquet(eqtl_inpt)
    qtl['chrom'] = chrom
    qtl['snp'] = qtl['variant_id'].apply(extract_snp)
    qtl['type'] = "quant"
    qtl['maf'] = qtl['af'].apply(lambda af: 1 - af if af > 0.5 else af)
    qtl['n'] = (qtl['ma_count'] / (2 * qtl['maf'])).round()

    # Read in the GWAS summary statistics
    gwas_result = pyreadr.read_r(gwas_input)
    gwas = list(gwas_result.values())[0]

    print(f"Done reading/preparing data: {chrom}")

    results = []
    # Detect number of cores
    n_cores = 4
    print(f"Using {n_cores} cores for parallel processing.")
    # Process each candidate gene in parallel
    with ProcessPoolExecutor(max_workers=n_cores) as executor:
        # Submit all tasks
        futures = [
            executor.submit(colocWindow_for_gene, gene, qtl, gwas, N_eqtl, N_gwas, chrom, output_pre)
            for gene in sig_genes
        ]
        # Collect results with a progress bar
        for future in tqdm(as_completed(futures), total=len(futures), desc="Processing candidate genes"):
            res = future.result()
            if res is not None:
                results.append(res)
    
    # Combine results into a DataFrame
    fin_coloc = pd.DataFrame(results)
    print(f"Finish: {chrom}")
    return fin_coloc

def main():
    # Extract input
    args = parse_args()
    
    # Read the candidate genes from the shortList file
    sig_genes_df = pd.read_csv(args.shortList, sep="\t", compression='infer')
    sig_genes = sig_genes_df["phenotype_id"].tolist()
    
    dt = coloc_chrom(eqtl_inpt=args.eqtl,
                     gwas_input=args.gwas,
                     sig_genes=sig_genes,
                     chrom=args.chromosome,
                     N_eqtl=args.N_eqtl,
                     N_gwas=args.N_gwas,
                     output_pre=args.prefix)
    
    # Write results to output file
    output_file = f"coloc_eqtl_{args.prefix}_{args.chromosome}.txt"
    dt.to_csv(output_file, sep="\t", index=False)
    print(f"Output written to {output_file}")

if __name__ == '__main__':
    main()

