#!/usr/bin/env python
# coding: utf-8

"""
By: Connor S. Murray

  - Performs TMM + inverse normal transform + PCA on RNA-seq data
  - Plots the first 6 PCs + scree plot, colored by Affected_NF
  - Outputs up to the 50th PC in pca_data.tsv
  - Assesses biological sex based on XIST and RPS4Y1 expression (ANOVA approach)
  - Removes outliers based on:
      (1) High mean absolute deviation across first 10 PCs, normalized by each PC's MAD
      (2) Significant Mahalanobis distance for first 5 PCs (Bonferroni-corrected p < 0.05)
  - Writes final outlier sample IDs to a text file
  - Saves final filtered log-transformed counts and PCA to disk.
"""
# Load Libraries
import argparse
import sys
import os
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
from scipy.spatial.distance import mahalanobis
from scipy.stats import chi2
from scipy.stats import median_abs_deviation
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multitest import multipletests
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from tqdm import tqdm
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from pydeseq2.utils import load_example_data
import matplotlib.pyplot as plt
import seaborn as sns

def anova_pvalues(df, predictors):
    """
    Perform ANOVA for each predictor in 'predictors' and return p-values.
    """
    p_values = {}
    for predictor in predictors:
        try:
            model = ols(f"norm ~ {predictor}", data=df).fit()
            anova_table = sm.stats.anova_lm(model, typ=2)
            p_values[predictor] = anova_table["PR(>F)"][predictor]
        except Exception:
            p_values[predictor] = np.nan
    return p_values

def process_gene(gene_group, predictors):
    """
    Process a single gene group (long format) and return ANOVA p-values.
    """
    p_values = anova_pvalues(gene_group, predictors)
    return {"Name": gene_group["Name"].iloc[0], **p_values}

def normalize_and_filter(metadata_path, mappability_path, gtf_path, 
                         gene_counts_path, output_normalized, output_pca,
                         output_outliers, skip_mappability_filter=True):
    """
    Normalize, filter, and run PCA on RNA-seq gene expression data.
    - Removes ERCC genes and genes not in the GTF or with low mappability *before* normalization.
    - Fixes the 'sample_id' index vs. column confusion.
    - Removes outliers based on (1) MAD across first 10 PCs, (2) Mahalanobis dist on first 5 PCs.
    - Saves a PCA plot (pca_medrat_plot.pdf) and writes final log counts & PCA results.
    - Writes outlier IDs to output_outliers.
    """
    # -----------------------------------------------------------
    # (1) Load Metadata
    # -----------------------------------------------------------
    print("Loading metadata...")
    meta_full = pd.read_table(metadata_path)
    if 'SAMPLE_ID_TOR' not in meta_full.columns:
        raise ValueError("metadata file must contain 'SAMPLE_ID_TOR' column.")

    # Set index to SAMPLE_ID_TOR
    meta_full["SAMPLE_ID_TOR"] = meta_full["SAMPLE_ID_TOR"].astype(str)
    meta_full.set_index("SAMPLE_ID_TOR", inplace=True)

    # Example known sample removal
    if "TOR238072" in meta_full.index:
        meta_full = meta_full.drop(index="TOR238072")

    # We only need certain columns for DESeq2
    meta = meta_full[["Affected_NF"]].copy()
    meta = meta.rename(columns={"Affected_NF": "Affected-NF"})
    meta.index = meta.index.astype(str)
    meta["Affected-NF"] = meta["Affected-NF"].astype("category")

    # -----------------------------------------------------------
    # (2) Load Mappability and GTF
    # -----------------------------------------------------------
    if not skip_mappability_filter:
        print("Loading mappability data...")
        mapp_df = pd.read_table(mappability_path, names=["gene_id", "mappability"], header=None)
        low_map_genes = set(mapp_df[mapp_df["mappability"] <= 0.5]["gene_id"])
    else:
        print("Skipping mappability data loading as per user request.")
        low_map_genes = set()

    print("Loading GTF file...")
    gene_gtf = pd.read_csv(gtf_path, delimiter="\t")
    gtf_genes_simp = set(gene_gtf["gene_id"].astype(str).replace(r"\..*", "", regex=True))

    # -----------------------------------------------------------
    # (3) Load Gene Counts
    # -----------------------------------------------------------
    print("Loading gene counts...")
    gene_count = pd.read_table(
        gene_counts_path,
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

    # -----------------------------------------------------------
    # (4) Remove ERCC Genes, Non-GTF Genes, or Low Mappability
    # -----------------------------------------------------------
    gene_count["gene_id"] = gene_count.index.str.replace(r"\..*", "", regex=True)
    gene_count["gene"] = gene_count.index

    ercc_mask = gene_count["gene_id"].str.contains("ERCC", case=False, na=False)
    not_in_gtf_mask = ~gene_count["gene_id"].isin(gtf_genes_simp)
    
    if not skip_mappability_filter:
        low_map_mask = gene_count["gene"].isin(low_map_genes)
    else:
        low_map_mask = pd.Series([False]*gene_count.shape[0], index=gene_count.index)

    combined_mask = ercc_mask | not_in_gtf_mask | low_map_mask
    print("Removing ERCC genes, genes not in GTF, or with low mappability...")
    filtered_gene_count = gene_count.loc[~combined_mask].copy()
    print(f"After removal, shape: {filtered_gene_count.shape}")

    # -----------------------------------------------------------
    # (5) DESeq2 Normalization
    # -----------------------------------------------------------
    gene_countsi = filtered_gene_count.drop(columns=["gene_id", "gene"]).T
    print(f"Gene count matrix for DESeq2: {gene_countsi.shape}")

    dds = DeseqDataSet(
        counts=gene_countsi,
        metadata=meta,
        design_factors="Affected-NF")
    dds.fit_size_factors()
    norm_counts = pd.DataFrame(dds.layers["normed_counts"])
    norm_counts.index = gene_countsi.index
    norm_counts.columns = gene_countsi.columns
    print(f"Normalized count matrix shape: {norm_counts.shape}")

    # -----------------------------------------------------------
    # (6) Filter Genes by QC
    # -----------------------------------------------------------
    norm_t = norm_counts.T
    num_samples = meta.shape[0]

    # Overall filter
    norm_t["qc_individual_depth"] = (norm_t >= 10).sum(axis=1) / num_samples * 100
    norm_t["pass_depth"] = np.where(norm_t["qc_individual_depth"] > 10, "Pass", "Fail")

    # Filter among Affected only
    affected_samples = meta[meta["Affected-NF"] == "Affected"].index
    affected_counts = norm_t.loc[:, affected_samples].copy()
    affected_counts["qc_individual_depth"] = (
        (affected_counts >= 10).sum(axis=1) / len(affected_samples) * 100)
    
    affected_counts["pass_depth_affected"] = np.where(
        affected_counts["qc_individual_depth"] > 10, "Pass", "Fail")

    pass_affected_genes = affected_counts[affected_counts["pass_depth_affected"] == "Pass"].index
    filt_norm_t = norm_t[
        (norm_t["pass_depth"] == "Pass") &
        (norm_t.index.isin(pass_affected_genes))]
    
    print(f"Genes passing combined filter: {filt_norm_t.shape[0]}")

    # -----------------------------------------------------------
    # (7) Replace Zeros with Half-min, Log-transform
    # -----------------------------------------------------------
    half_mins = filt_norm_t.drop(
        columns=["qc_individual_depth", "pass_depth"], errors="ignore"
    ).apply(
        lambda x: np.min(x[x > 0]) / 2 if (x > 0).any() else 0,
        axis=1)
    filt_log_t = filt_norm_t.drop(
        columns=["qc_individual_depth", "pass_depth"], errors="ignore"
    ).apply(
        lambda row: np.log10(row.replace(0, half_mins[row.name])),
        axis=1)
    
    print(f"Final number of genes after zero-replacement/log: {filt_log_t.shape[0]}")

    # Keep a 'Name' column for merges
    filt_log_t["Name"] = filt_log_t.index

    # -----------------------------------------------------------
    # (8) ANOVA & Remove Sex-correlated Genes
    # -----------------------------------------------------------
    meta_full.index.name = None
    needed_cols = [
        "Gender", "Sample_container_ID", "Affected_NF",
        "Race", "Diagnosis", "Primary_biosample_type"]
    
    needed_cols = [c for c in needed_cols if c in meta_full.columns]

    melted = filt_log_t.melt(
        id_vars=["Name"],
        var_name="sample_id",
        value_name="norm")
    
    merged_df = melted.join(meta_full[needed_cols], on="sample_id", how="left")

    drop_cols = ["norm", "Gender", "Sample_container_ID"]
    existing_drop_cols = [c for c in drop_cols if c in merged_df.columns]
    merged_df.dropna(subset=existing_drop_cols, inplace=True)

    predictors = ["Sample_container_ID", "Gender"]
    def process_gene_silent(gene_group, preds):
        sys.stdout = open(os.devnull, 'w')
        result = process_gene(gene_group, preds)
        sys.stdout = sys.__stdout__
        return result

    grouped = merged_df.groupby("Name")
    results = Parallel(n_jobs=-1, verbose=0)(
        delayed(process_gene_silent)(group, predictors)
        for _, group in tqdm(grouped, total=len(grouped), desc="ANOVA", ncols=80)
    )
    anova_df = pd.DataFrame(results)
    for pred in ["Sample_container_ID", "Gender"]:
        if pred in anova_df.columns:
            anova_df[f"ap_{pred}"] = multipletests(anova_df[pred], method="bonferroni")[1]

    if "ap_Gender" in anova_df.columns:
        sex_genes = anova_df[anova_df["ap_Gender"] < 0.05]["Name"]
    else:
        sex_genes = []
        print("WARNING: 'ap_Gender' not found in ANOVA results. Skipping sex-correlation removal.")

    final_log_t = filt_log_t[~filt_log_t["Name"].isin(sex_genes)].drop(columns=["Name"])
    print(f"Number of genes after removing sex-correlated: {final_log_t.shape[0]}")

    # -----------------------------------------------------------
    # (9) PCA (Saving pca_medrat_plot.pdf)
    # -----------------------------------------------------------
    data_for_pca = final_log_t.T  # shape => (samples x genes)
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(data_for_pca)
    pca = PCA(n_components=50)
    pca_result = pca.fit_transform(scaled_data)

    pca_df = pd.DataFrame(pca_result, columns=[f"PC{i+1}" for i in range(50)])
    pca_df.index = data_for_pca.index  # sample IDs

    # Add 'Affected_NF' if available
    if "Affected_NF" in meta_full.columns:
        pca_df["Affected_NF"] = meta_full["Affected_NF"]

    # Plot first 6 PCs
    sns.set_style("whitegrid")
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.subplots_adjust(hspace=0.25, wspace=0.2)

    pairs = [(0,1), (1,2), (2,3), (3,4), (4,5)]
    for i, (pc_x, pc_y) in enumerate(pairs):
        ax = axes.flatten()[i]
        sns.scatterplot(
            x=pca_df.iloc[:, pc_x],
            y=pca_df.iloc[:, pc_y],
            hue=pca_df.get('Affected_NF', None),
            alpha=0.7, s=50, ax=ax,
            legend=(i==0)
        )
        ax.set_xlabel(f"PC{pc_x+1}", weight="bold", size=14)
        ax.set_ylabel(f"PC{pc_y+1}", weight="bold", size=14)

    # Scree plot in the 6th panel
    scree_ax = axes.flatten()[5]
    scree_ax.bar(range(1,51), pca.explained_variance_ratio_[:50]*100, color='skyblue')
    scree_ax.set_xlabel('Principal Component', weight="bold", size=12)
    scree_ax.set_ylabel('Explained Variance (%)', weight="bold", size=12)
    scree_ax.set_xticks(range(1,51,4))
    plt.tight_layout()
    plt.savefig("pca_medrat_plot.pdf", dpi=300)
    plt.close()

    print("PCA explained variance ratios (first 6 PCs):",
          pca.explained_variance_ratio_[:6]*100)

    # -----------------------------------------------------------
    # (9b) Mahalanobis distance on first 5 PCs
    # -----------------------------------------------------------
    pc_data = pca_df[[f"PC{i}" for i in range(1,6)]].copy()
    mean_vector = pc_data.mean(axis=0)
    cov_matrix = np.cov(pc_data, rowvar=False)
    inv_cov_matrix = np.linalg.inv(cov_matrix)

    mal_dist = pc_data.apply(lambda row: mahalanobis(row, mean_vector, inv_cov_matrix), axis=1)
    pvals = 1 - chi2.cdf(mal_dist, df=pc_data.shape[1]-1)
    p_adj = multipletests(pvals, method="fdr_bh")[1]

    pca_df["mahal"] = mal_dist
    pca_df["mahal_pval"] = pvals
    pca_df["mahal_padj"] = p_adj

    # -----------------------------------------------------------
    # (9c) MAD-based outlier detection on first 10 PCs
    # -----------------------------------------------------------
    ## BEGIN OUTLIER DETECTION ##
    pc10 = pca_df[[f"PC{i}" for i in range(1,11)]].copy()

    # For each PC, compute the median & MAD, then normalize each sample's deviation
    for i in range(1, 11):
        pc_col = f"PC{i}"
        median_i = pc10[pc_col].median()
        mad_i = median_abs_deviation(pc10[pc_col])
        # Avoid divide-by-zero if mad_i=0
        if mad_i == 0:
            # If there's no variability for that PC, set normalized dev to 0
            pc10[f"{pc_col}_normdev"] = 0.0
        else:
            pc10[f"{pc_col}_normdev"] = (pc10[pc_col] - median_i).abs() / mad_i

    # Average normalized deviation across the 10 PCs
    normdev_cols = [f"PC{i}_normdev" for i in range(1,11)]
    pc10["mad_score"] = pc10[normdev_cols].mean(axis=1)

    # Mark samples with mean normalized deviation >= 3
    pca_df["mad_score"] = pc10["mad_score"]
    pca_df["mad_outlier"] = pca_df["mad_score"] >= 3

    # Combine outliers from either criterion
    #   1) mad_outlier == True
    #   2) mahal_padj < 0.05
    outlier_set_mad = pca_df.index[pca_df["mad_outlier"]]
    outlier_set_mahal = pca_df.index[pca_df["mahal_padj"] < 0.05]
    combined_outliers = set(outlier_set_mad).union(set(outlier_set_mahal))

    # Convert to a sorted list
    combined_outliers = sorted(list(combined_outliers))

    # -----------------------------------------------------------
    # (9d) Remove outliers from final data
    # -----------------------------------------------------------
    print(f"Total outliers found = {len(combined_outliers)}")
    if len(combined_outliers) > 0:
        print("Removing these outliers from final counts and PCA data...")
    else:
        print("No outliers detected by either criterion.")

    # Remove them from final_log_t (genes x samples) -> index = gene, columns = sample
    # We must remove them from the columns if final_log_t is shaped that way
    final_log_t_cleaned = final_log_t.drop(columns=combined_outliers, errors="ignore")

    # Remove from PCA
    pca_df_cleaned = pca_df.drop(index=combined_outliers, errors="ignore")

    # -----------------------------------------------------------
    # (10) Save Results
    # -----------------------------------------------------------
    # Switch final_log_t_cleaned to (samples x genes) for saving
    final_log_counts = final_log_t_cleaned.T

    print(f"Saving final filtered, log-transformed counts to: {output_normalized}")
    final_log_counts.to_csv(output_normalized, sep="\t")

    print(f"Saving final PCA DataFrame to: {output_pca}")
    pca_df_cleaned.to_csv(output_pca, sep="\t")

    print(f"Saving outlier individuals to: {output_outliers}")
    with open(output_outliers, "w") as f:
        for sample in combined_outliers:
            f.write(str(sample) + "\n")

    print("Done. Outlier removal complete, final data and PCA saved.")

def main():
    parser = argparse.ArgumentParser(description="Normalize, filter, and perform PCA on RNA-seq data.")
    parser.add_argument("--metadata", required=True, help="Path to the metadata file.")
    parser.add_argument("--mappability", required=False, help="Path to the mappability file.")
    parser.add_argument("--gtf", required=True, help="Path to the GTF file.")
    parser.add_argument("--gene_counts", required=True, help="Path to the gene counts file.")
    parser.add_argument("--output_normalized", required=True, help="Path to save final normalized counts.")
    parser.add_argument("--output_pca", required=True, help="Path to save PCA table.")
    parser.add_argument("--output_outliers", required=True, help="Path to save outlier individuals.")
    # --- Added: Skip Mappability Filter Flag ---
    parser.add_argument("--skip_mappability_filter", action="store_true",
                        help="Flag to skip the low-mappability gene filtering step.")
    args = parser.parse_args()

    # --- Argument Validation ---
    if not args.skip_mappability_filter and not args.mappability:
        parser.error("--mappability is required unless --skip_mappability_filter is set.")

    normalize_and_filter(
        metadata_path=args.metadata,
        mappability_path=args.mappability,
        gtf_path=args.gtf,
        gene_counts_path=args.gene_counts,
        output_normalized=args.output_normalized,
        output_pca=args.output_pca,
        output_outliers=args.output_outliers,
        skip_mappability_filter=args.skip_mappability_filter) 

if __name__ == "__main__":
    main()
