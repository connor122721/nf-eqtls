#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Process for analyzing eqtl output
process analysis_eqtl_saturation {

    // Publish the output to the specified directory
    shell = '/usr/bin/env bash'
    publishDir "${params.out}/analyses", mode: 'copy'

    input:
        path metadata
        path gtf

    output:
        val "${best_k}"

    script:
        """
        module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1

        Rscript ${params.scripts_dir}/analysis_eQTL_saturation.R \\
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

        echo "Finish"
        """
}