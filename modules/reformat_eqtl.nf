#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//  Process for reformatting eqtl

// Process for reformatting eqtl
process reformat_eqtl {

    // Publish the output to the specified directory
    shell = '/usr/bin/env bash'
    publishDir "${params.out}/eqtl", mode: 'copy'

    input:
        path metadata
        path gtf
        path related_individuals
        path rna_outliers
        path dna_outliers
        path norm_gene_count
        path pca_tmm
        path tensorqtl_pca

    output:
        path("topchef_cov_RNApc1_*.txt")
        path("filt_rnaseq_norm_topchef_unrelated.bed")
        path("filt_rnaseq_norm_topchef_unrelated.sort.bed")
        path("topchef_samples_1_15_25.txt"), emit: sample_list
        path("rna_pca_elbow_best_k")

    script:
        """
        module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1

        Rscript ${params.scripts_dir}/reformat_eqtl.R \\
            --metadata ${metadata} \\
            --gene_gtf ${gtf} \\
            --related_individuals ${related_individuals} \\
            --rna_outliers ${rna_outliers} \\
            --dna_outliers ${dna_outliers} \\
            --norm_tmm ${norm_gene_count} \\
            --pca_tmm ${pca_tmm} \\
            --pca_snp ${tensorqtl_pca}

        # Sort bed file
        sort -k1,1 -k2,2n filt_rnaseq_norm_topchef_unrelated.bed > \\
            filt_rnaseq_norm_topchef_unrelated.sort.bed
        """
}