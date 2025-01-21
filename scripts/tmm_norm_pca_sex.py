#!/usr/bin/env python
# coding: utf-8

"""
By: Connor S. Murray
TMM + inverse normal transform + PCA script
Plot the first 6 PCs + scree plot, color by Affected_NF
and output up to the 30th PC in pca_data.tsv
Additionally, assess biological sex based on XIST and RPS4Y1 expression.
"""

# Load Libraries
import argparse
import sys
import os
import pandas as pd
import numpy as np
from scipy.stats import rankdata, norm
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
import qtl.norm  # Must have this installed to use qtl.norm.edger_cpm

def parse_arguments():
    parser = argparse.ArgumentParser(description="TMM + INT normalization + PCA pipeline (up to PC30) with 6-PC plot and sex assessment.")
    parser.add_argument("--metadata", required=False, help="Path to metadata file (must include SAMPLE_ID_TOR, Affected_NF).")
    parser.add_argument("--gene_counts", required=False, help="Path to gene counts file (GCT format).")
    parser.add_argument("--gtf", required=False, help="Path to GTF bed file (e.g., gencode.v34...).")
    parser.add_argument("--mappability", required=False, help="Path to gzipped mappability file (hg38_gene_mappability.txt.gz).")
    parser.add_argument("--output_norm", default="normalized_data_tmm.tsv",
                        help="Output path for final TMM + inverse-normal-transformed data (genes x samples).")
    parser.add_argument("--output_pca", default="pca_data_tmm.tsv",
                        help="Output path for PCA data (PC1..PC30).")
    parser.add_argument("--output_plot_pdf", default="pca_plot_tmm.pdf",
                        help="Output path for PCA plot PDF (first 6 PCs + scree).")
    parser.add_argument("--output_sex_plot_pdf", default="sex_assessment_plot.pdf",
                        help="Output path for Sex Assessment plot PDF.")
    parser.add_argument("--low_expression_threshold", type=float, default=0.1,
                        help="Minimum CPM threshold to consider a gene expressed.")
    parser.add_argument("--sample_expression_frac", type=float, default=0.3,
                        help="Fraction of samples that must exceed the CPM threshold for a gene to be retained.")
    parser.add_argument("--affected_expression_frac", type=float, default=0.3,
                        help="Fraction of affected samples that must exceed the coverage threshold for a gene to be retained.")
    parser.add_argument("--xist_gene", default="ENSG00000229807",
                        help="Ensembl Gene ID for XIST (default: ENSG00000229807).")
    parser.add_argument("--rps4y1_gene", default="ENSG00000129824",
                        help="Ensembl Gene ID for RPS4Y1 (default: ENSG00000129824).")
    parser.add_argument("--xist_threshold", type=float, default=10.0,
                        help="CPM threshold to consider XIST as highly expressed (default: 10).")
    parser.add_argument("--rps4y1_threshold", type=float, default=1.0,
                        help="CPM threshold to consider RPS4Y1 as highly expressed (default: 1).")
    parser.add_argument("--exclude_unclear_sex", action="store_true",
                        help="Flag to exclude samples with unclear sex from downstream analyses.")
    parser.add_argument("--test", action="store_true",
                        help="Run the script in testing mode with predefined file paths.")
    return parser.parse_args()

def load_data(gene_counts_path, metadata_path):
    print("Loading TPM-normalized gene counts...")
    gene_counts = pd.read_table(gene_counts_path, compression="gzip", skiprows=2, index_col="Name")
    print(f"Gene counts shape: {gene_counts.shape}")

    print("Loading metadata...")
    meta_full = pd.read_table(metadata_path)
    if "SAMPLE_ID_TOR" not in meta_full.columns:
        raise ValueError("Metadata file must contain 'SAMPLE_ID_TOR' column.")

    # Index by SAMPLE_ID_TOR
    meta_full["SAMPLE_ID_TOR"] = meta_full["SAMPLE_ID_TOR"].astype(str)
    meta_full.set_index("SAMPLE_ID_TOR", inplace=True)

    # Drop sample if needed
    if "TOR238072" in meta_full.index:
        meta_full.drop("TOR238072", inplace=True)

    # We expect 'Affected_NF' for coloring
    if "Affected_NF" not in meta_full.columns:
        raise ValueError("Metadata must have 'Affected_NF' column for coloring PCA points.")

    # Minimal usage for alignment
    meta = meta_full[["Affected_NF"]].copy()
    meta.rename(columns={"Affected_NF": "Affected-NF"}, inplace=True)
    meta["Affected-NF"] = meta["Affected-NF"].astype("category")

    # Align columns with metadata
    gene_counts.columns = gene_counts.columns.astype(str)
    common_samples = meta.index.intersection(gene_counts.columns)
    meta = meta.loc[common_samples]
    gene_counts = gene_counts.loc[:, common_samples]
    print(f"Number of samples after alignment: {len(common_samples)}")

    return gene_counts, meta, meta_full

# Function to infer sex based on XIST and RPS4Y1 normalized expression
def infer_sex(ccm, meta, xist_gene, rps4y1_gene, xist_threshold, rps4y1_threshold):
    print("Inferring sex based on XIST and RPS4Y1 expression...")

    # Handle possible version suffixes in gene IDs (e.g., ENSG00000229807.1)
    base_xist_gene = xist_gene.split('.')[0]
    base_rps4y1_gene = rps4y1_gene.split('.')[0]

    # Find all gene IDs that match the base IDs
    xist_matches = ccm.index[ccm.index.str.startswith(base_xist_gene)]
    rps4y1_matches = ccm.index[ccm.index.str.startswith(base_rps4y1_gene)]

    if len(xist_matches) == 0:
        raise ValueError(f"No matches found for XIST gene ID: {xist_gene}")
    if len(rps4y1_matches) == 0:
        raise ValueError(f"No matches found for RPS4Y1 gene ID: {rps4y1_gene}")

    # Aggregate CPM values for XIST and RPS4Y1 by summing across variants
    xist_cpm = ccm.loc[xist_matches].sum(axis=0)
    rps4y1_cpm = ccm.loc[rps4y1_matches].sum(axis=0)

    # Initialize sex column
    inferred_sex = []

    for sample in ccm.columns:
        xist = xist_cpm[sample]
        rps4y1 = rps4y1_cpm[sample]

        if (xist >= xist_threshold) and (rps4y1 <= rps4y1_threshold):
            inferred_sex.append("Female")
        elif (xist <= xist_threshold) and (rps4y1 >= rps4y1_threshold):
            inferred_sex.append("Male")
        else:
            inferred_sex.append("Unclear")

    meta = meta.copy()
    meta["Inferred_Sex"] = inferred_sex
    num_females = (meta["Inferred_Sex"] == "Female").sum()
    num_males = (meta["Inferred_Sex"] == "Male").sum()
    num_unclear = (meta["Inferred_Sex"] == "Unclear").sum()
    print(f"Inferred Sex Counts: Female={num_females}, Male={num_males}, Unclear={num_unclear}")

    return meta

# Function to plot biological sex assessment
def plot_sex_assessment(ccm, meta, xist_gene, rps4y1_gene, output_plot_pdf, xist_threshold, rps4y1_threshold):
    print("Generating sex assessment plot...")
    # Handle possible version suffixes
    base_xist_gene = xist_gene.split('.')[0]
    base_rps4y1_gene = rps4y1_gene.split('.')[0]

    # Find all gene IDs that match the base IDs
    xist_matches = ccm.index[ccm.index.str.startswith(base_xist_gene)]
    rps4y1_matches = ccm.index[ccm.index.str.startswith(base_rps4y1_gene)]

    # Aggregate CPM values by summing across variants
    xist_cpm = ccm.loc[xist_matches].sum(axis=0)
    rps4y1_cpm = ccm.loc[rps4y1_matches].sum(axis=0)

    plot_df = pd.DataFrame({
        xist_gene: xist_cpm,
        rps4y1_gene: rps4y1_cpm,
        "Inferred_Sex": meta["Inferred_Sex"]})

    plt.figure(figsize=(8, 8))
    sns.scatterplot(
        data=plot_df,
        x=xist_gene,
        y=rps4y1_gene,
        hue="Inferred_Sex",
        palette={"Female": "red", "Male": "blue", "Unclear": "grey"},
        alpha=0.7,
        s=100)
    plt.axvline(x=xist_threshold, color='grey', linestyle='--', linewidth=1)
    plt.axhline(y=rps4y1_threshold, color='grey', linestyle='--', linewidth=1)
    plt.title("Sex Assessment Based on XIST and RPS4Y1 Expression", fontsize=16)
    plt.xlabel(f"{xist_gene} CPM", fontsize=14)
    plt.ylabel(f"{rps4y1_gene} CPM", fontsize=14)
    plt.legend(title="Inferred Sex")
    plt.tight_layout()

    plt.savefig(output_plot_pdf, dpi=300)
    plt.close()
    print(f"Sex assessment plot saved as {output_plot_pdf}")

# Function to filter samples based on biological sex assessment
def filter_samples_based_on_sex(ccm, meta, exclude_unclear_sex):
    if exclude_unclear_sex:
        print("Excluding samples with unclear sex...")
        initial_sample_count = meta.shape[0]
        meta_filtered = meta[meta["Inferred_Sex"] != "Unclear"]
        ccm_filtered = ccm.loc[:, meta_filtered.index]
        excluded_samples = initial_sample_count - meta_filtered.shape[0]
        print(f"Number of samples excluded: {excluded_samples}")
        return ccm_filtered, meta_filtered
    else:
        print("Keeping all samples regardless of sex assessment.")
        return ccm, meta

def main():
    args = parse_arguments()

    # -------------------------------------------------------------------
    # 1) Load gene counts and metadata
    # -------------------------------------------------------------------
    gene_counts, meta, meta_full = load_data(args.gene_counts, args.metadata)

    # -------------------------------------------------------------------
    # 2) TMM normalization
    # -------------------------------------------------------------------
    print("Performing TMM normalization (edger_cpm)...")
    ccm = qtl.norm.edger_cpm(gene_counts)
    print("Example of normalized CPM values:")
    print(ccm.head())
    print(f"Pre-filtering genes: {ccm.shape}")

    # -------------------------------------------------------------------
    # 3) Assess Biological Sex Based on XIST and RPS4Y1 Expression
    # -------------------------------------------------------------------
    meta = infer_sex(ccm=ccm,
                     meta=meta,
                     xist_gene=args.xist_gene,
                     rps4y1_gene=args.rps4y1_gene,
                     xist_threshold=args.xist_threshold,
                     rps4y1_threshold=args.rps4y1_threshold)

    # Generate Sex Assessment Plot
    plot_sex_assessment(ccm=ccm,
                        meta=meta,
                        xist_gene=args.xist_gene,
                        rps4y1_gene=args.rps4y1_gene,
                        output_plot_pdf=args.output_sex_plot_pdf,
                        xist_threshold=args.xist_threshold,
                        rps4y1_threshold=args.rps4y1_threshold)

    # Optionally, exclude samples with unclear sex
    ccm, meta = filter_samples_based_on_sex(ccm, meta, args.exclude_unclear_sex)

    # -------------------------------------------------------------------
    # 4) Filter lowly expressed genes
    # -------------------------------------------------------------------
    print("Filtering lowly expressed genes...")
    low_expr_thr = args.low_expression_threshold
    frac_samples = args.sample_expression_frac
    sample_threshold = int(np.ceil(frac_samples * ccm.shape[1]))
    gene_filter = (ccm > low_expr_thr).sum(axis=1) >= sample_threshold
    ccm_filtered = ccm[gene_filter]
    print(f"Genes retained after expression filter: {ccm_filtered.shape[0]}")

    # -------------------------------------------------------------------
    # 5) Load GTF and remove non-GTF genes
    # -------------------------------------------------------------------
    print("Loading GTF file...")
    gtf_df = pd.read_csv(args.gtf, delimiter="\t", header=None,
                         names=["chrom", "start", "end", "gene_id"])
    gtf_ids = set(gtf_df["gene_id"].astype(str))
    gene_in_gtf = ccm_filtered.index.isin(gtf_ids)
    ccm_in_gtf = ccm_filtered[gene_in_gtf]
    print(f"Genes after GTF filter: {ccm_in_gtf.shape[0]} (out of {ccm_filtered.shape[0]})")

    # -------------------------------------------------------------------
    # 6) Remove ERCC genes
    # -------------------------------------------------------------------
    ercc_mask = ccm_in_gtf.index.str.contains("ERCC", case=False, na=False)
    ercc_count = ercc_mask.sum()
    print(f"Number of ERCC genes: {ercc_count}")
    ccm_in_gtf = ccm_in_gtf[~ercc_mask]
    print(f"Shape after removing ERCC: {ccm_in_gtf.shape}")

    # -------------------------------------------------------------------
    # 7) Remove low-mappability genes
    # -------------------------------------------------------------------
    print("Loading mappability data...")
    mapp_df = pd.read_table(args.mappability, header=0, names=["gene_id", "mappability"])
    low_map = set(mapp_df.loc[mapp_df["mappability"] <= 0.5, "gene_id"])
    print(f"Number of low-mappability genes: {len(low_map)}")
    keep_mask = ~ccm_in_gtf.index.isin(low_map)
    ccm_final = ccm_in_gtf[keep_mask]
    print(f"Final shape after removing low-mappability genes: {ccm_final.shape}")

    # -------------------------------------------------------------------
    # 8) Filter genes with low coverage among affected samples
    # -------------------------------------------------------------------
    print("Filtering genes based on coverage in affected samples...")
    # Identify affected samples
    affected_samples = meta[meta["Affected-NF"] == "Affected"].index
    if len(affected_samples) == 0:
        raise ValueError("No samples labeled as 'Affected' found in metadata.")

    # Extract counts for affected samples
    affected_counts = ccm_final.loc[:, affected_samples].copy()

    # Calculate the percentage of affected samples with coverage >= threshold per gene
    coverage_frac_threshold = args.affected_expression_frac
    genes_meeting_coverage = (affected_counts > low_expr_thr).sum(axis=1) / len(affected_samples) >= coverage_frac_threshold

    # Apply the coverage filter
    ccm_filtered_final = ccm_final[genes_meeting_coverage]
    print(f"Genes retained after coverage filter: {ccm_filtered_final.shape[0]}")

    # -------------------------------------------------------------------
    # 9) Inverse normal transform (per gene)
    # -------------------------------------------------------------------
    print("Applying inverse normal transform (per gene)...")

    def inverse_normal_transform(row):
        # row => 1D array of length = number of samples
        ranks = rankdata(row, method='average')
        quantiles = ranks / (len(ranks) + 1.0)
        return pd.Series(norm.ppf(quantiles), index=row.index)

    ccm_tmm_int = ccm_filtered_final.apply(inverse_normal_transform, axis=1)
    print(f"Final shape after INT transform: {ccm_tmm_int.shape}")
    if ccm_tmm_int.isnull().sum().sum() > 0:
        raise ValueError("Transformed data contains NaN values after INT step!")

    # -------------------------------------------------------------------
    # 10) Principal component analysis
    # -------------------------------------------------------------------
    # shape => (samples x genes)
    data_for_pca = ccm_tmm_int.T.astype(float)
    print(f"Performing PCA with 30 components on shape: {data_for_pca.shape}")
    pca = PCA(n_components=30)
    data_scaled = StandardScaler().fit_transform(data_for_pca)
    pca_result = pca.fit_transform(data_scaled)

    # Build a DataFrame
    pc_columns = [f"PC{i+1}" for i in range(30)]
    pca_df = pd.DataFrame(pca_result, columns=pc_columns, index=data_for_pca.index)

    # Add Affected_NF for coloring
    # meta has index = sample IDs => reindex to align with pca_df
    pca_df["Affected_NF"] = meta.reindex(pca_df.index)["Affected-NF"]
    print("PCA shape:", pca_df.shape)

    # -------------------------------------------------------------------
    # 11) Plot the first 6 PCs + scree
    # -------------------------------------------------------------------
    print(f"Plotting first 6 PCs and scree to {args.output_plot_pdf}...")
    sns.set_style("whitegrid")
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.subplots_adjust(hspace=0.25, wspace=0.2)

    # Scatter plot pairs: PC1 vs PC2, etc.
    pairs = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)]

    for i, (pcx, pcy) in enumerate(pairs):
        ax = axes.flat[i]
        sns.scatterplot(
            x=pca_df[pc_columns[pcx]],
            y=pca_df[pc_columns[pcy]],
            hue=pca_df["Affected_NF"],
            alpha=0.7, s=50,
            ax=ax, legend=(i == 0)  # show legend only on first scatter
        )
        ax.set_xlabel(pc_columns[pcx], fontweight="bold", fontsize=12)
        ax.set_ylabel(pc_columns[pcy], fontweight="bold", fontsize=12)

    # Scree plot in the last subplot (axes.flat[5])
    scree_ax = axes.flat[5]
    explained_var = pca.explained_variance_ratio_ * 100  # convert to percentage
    scree_ax.bar(range(1, 31), explained_var[:30], color="skyblue")
    scree_ax.set_xlabel("Principal Component", fontweight="bold", fontsize=12)
    scree_ax.set_ylabel("Explained Variance (%)", fontweight="bold", fontsize=12)
    scree_ax.set_xticks(range(1, 31, 4))  # label every 4th PC
    scree_ax.set_ylim(0, max(explained_var[:30]) * 1.2)

    plt.tight_layout()
    plt.savefig(args.output_plot_pdf, dpi=300)
    plt.close()
    print(f"Saved PCA plot as {args.output_plot_pdf}")

    # -------------------------------------------------------------------
    # 12) Save Outputs
    # -------------------------------------------------------------------
    # ccm_tmm_int is (genes x samples), we'll keep that for normalized_data
    print(f"Saving final normalized data to: {args.output_norm}")
    ccm_tmm_int.to_csv(args.output_norm, sep="\t")

    # pca_df has PC1..PC30 + Affected_NF
    print(f"Saving PCA data (PC1..PC30) to: {args.output_pca}")
    pca_df.to_csv(args.output_pca, sep="\t")

    print("Done: TMM + INT normalization, sex assessment, 30-PC PCA, and plots.")

if __name__ == "__main__":
    main()
