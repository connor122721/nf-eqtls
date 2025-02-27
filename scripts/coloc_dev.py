#!/usr/bin/env python3

# Does coloc on given candidate gene list and processed GWAS
# Connor Murray; 2.26.2025

# Import packages
import argparse
import pandas as pd
import numpy as np
import re
import pyreadr  # for reading RDS files
import pyarrow.parquet as pq  # alternatively, you can use pd.read_parquet
from qtl import abf  # assuming a package “qtl” with coloc_abf exists

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

def coloc_chrom(eqtl_inpt, gwas_input, sig_genes, chrom, N_eqtl, N_gwas, output_pre):
    print(f"Running: {chrom}")
    
    # Read in the cis-eQTL data from a parquet file
    qtl = pd.read_parquet(eqtl_inpt)
    qtl['chrom'] = chrom
    qtl['snp'] = qtl['variant_id'].apply(extract_snp)
    qtl['type'] = "quant"
    qtl['maf'] = qtl['af'].apply(lambda af: 1 - af if af > 0.5 else af)
    qtl['n'] = (qtl['ma_count'] / (2 * qtl['maf'])).round()
    
    # Read in the GWAS summary statistics (RDS file)
    # pyreadr returns a dict; we assume the first object is our dataframe.
    gwas_result = pyreadr.read_r(gwas_input)
    gwas = list(gwas_result.values())[0]
    
    print(f"Done reading/preparing data: {chrom}")
    
    def colocWindow(dt1, focalGene, dt2, N1, N2):
        try:
            # Filter to the gene of interest and remove duplicate SNPs
            dt1_gene = dt1[dt1['phenotype_id'] == focalGene].drop_duplicates(subset=['variant_id'])
            if dt1_gene.empty:
                return None
            min_window = dt1_gene['snp'].min()
            max_window = dt1_gene['snp'].max()
            
            # Restrict dt2 to SNPs found in dt1_gene and vice versa
            dt2 = dt2[dt2['snpID_hg38'].isin(dt1_gene['variant_id'])]
            dt1_gene = dt1_gene[dt1_gene['variant_id'].isin(dt2['snpID_hg38'])]
            if dt1_gene.empty or dt2.empty:
                return None
            
            # Prepare dataset for coloc analysis (TOPChef/eQTL data)
            dataset1 = {
                "snp": dt1_gene["variant_id"].tolist(),
                "position": dt1_gene["snp"].tolist(),
                "type": "quant",
                "N": N1,
                "MAF": dt1_gene["maf"].tolist(),
                "pvalues": dt1_gene["pval_nominal"].tolist(),
                "beta": dt1_gene["slope"].tolist(),
                "varbeta": (dt1_gene["slope_se"] ** 2).tolist()
            }
            
            # Prepare dataset for coloc analysis (GWAS data)
            dataset2 = {
                "snp": dt2["snpID_hg38"].tolist(),
                "position": dt2["pos_hg38"].tolist(),
                "type": "quant",
                "N": N2,
                "MAF": dt2["maf"].tolist(),
                "pvalues": dt2["p_value"].tolist(),
                "beta": dt2["beta"].tolist(),
                "varbeta": (dt2["standard_error"] ** 2).tolist()
            }
            
            # Run colocalization analysis
            coloc_res = abf(dataset1=dataset1, dataset2=dataset2)
            
            # Extract the SNP with the maximum PP.H4
            results_df = coloc_res["results"]
            maxi = results_df["SNP.PP.H4"].max()
            snpi = results_df.loc[results_df["SNP.PP.H4"] == maxi, "snp"].iloc[0]
            
            # Collect results into a dictionary
            co = {
                "chrom": chrom,
                "minPos": min_window,
                "maxPos": max_window,
                "gene": focalGene,
                "min_p.eqtl": min(dataset1["pvalues"]),
                "min_p.gwas": min(dataset2["pvalues"]),
                "nsnps": coloc_res["summary"]["nsnps"],
                "PP.H0": coloc_res["summary"]["PP.H0.abf"],
                "PP.H1": coloc_res["summary"]["PP.H1.abf"],
                "PP.H2": coloc_res["summary"]["PP.H2.abf"],
                "PP.H3": coloc_res["summary"]["PP.H3.abf"],
                "PP.H4": coloc_res["summary"]["PP.H4.abf"],
                "H4_H3_ratio": coloc_res["summary"]["PP.H4.abf"] / coloc_res["summary"]["PP.H3.abf"],
                "maxSNP": snpi,
                "maxPP.H4": maxi,
                "gwas_pre": output_pre
            }
            return co
        except Exception as e:
            print(f"Error processing gene {focalGene}: {e}")
            return None
    
    # Iterate over each candidate gene and run colocalization
    results = []
    for gene in sig_genes:
        res = colocWindow(qtl, gene, gwas, N_eqtl, N_gwas)
        if res is not None:
            results.append(res)
    
    # Combine results into a DataFrame
    fin_coloc = pd.DataFrame(results)
    print(f"Finish: {chrom}")
    return fin_coloc

def main():
    args = parse_args()
    
    # Read the candidate genes from the shortList file (assuming tab-delimited, possibly gzipped)
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
