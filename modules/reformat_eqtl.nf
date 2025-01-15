#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// --------------------------------------------------
//  Processes related to covariate reformatting
// --------------------------------------------------

// Process for reformatting the covariates - split by group
process reformat_eqtl {

    // Publish the output to the specified directory
    shell = '/usr/bin/env bash'
    publishDir "${params.out}/eqtl", mode: 'copy'

    // Define the input and output
    input:
        path metadata
        path gene_gtf
        path related_individuals
        path samples_file
        path norm_tmm
        path pca_tmm

    output:
        path "topchef_cov_RNApc*.txt"
        path "filt_rnaseq_norm_topchef_unrelated.bed"
        path "topchef_samples*.txt"

    script:
        """
        module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1

        # Run the R script
        Rscript ${params.scripts_dir}/reformat_eqtl.R \\
            --metadata ${metadata} \\
            --gene_gtf ${gene_gtf} \\
            --related_individuals ${related_individuals} \\
            --norm_tmm ${norm_tmm} \\
            --pca_tmm ${pca_tmm}
        """
}
