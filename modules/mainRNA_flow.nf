#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// --------------------------------------------------
//  Processes related to RNA flow
// --------------------------------------------------

// Process for QC normalizing and PCA
process normalize_and_pca {
    
    shell = '/usr/bin/env bash'
    publishDir "${params.out}/rna", mode: 'copy'

    input:
        path meta_file
        path mapp_file
        path gtf_file
        path gene_count_file

    output:
        path "norm_medrat.tsv"
        path "pca_medrat.tsv"
        path "pca_medrat_plot.pdf"
        path "${params.rna_outliers}", emit: rna_outliers
        path "medratio.log"

    script:
        """
        module load miniforge/24.3.0-py3.11
        source activate base

        # Run script        
        python ${params.scripts_dir}/medratio_norm_pca.py \\
            --metadata ${meta_file} \\
            --gtf ${gtf_file} \\
            --gene_counts ${gene_count_file} \\
            --mappability ${mapp_file} \\
            --output_normalized norm_medrat.tsv \\
            --output_pca pca_medrat.tsv \\
            --output_outliers ${params.rna_outliers}
        # --skip_mappability_filter 

        cp .command.log medratio.log
        """
}

// Tmm normalization and PCA
process tmm_pipeline {
    
    shell = '/usr/bin/env bash'
    publishDir "${params.out}/rna", mode: 'copy'

    input:
        path meta_file
        path mapp_file
        path gtf_file
        path gene_count_file

    output:
        path "norm_tmm.tsv", emit: norm_gene_count
        path "pca_tmm.tsv", emit: pca_tmm
        path "pca_tmm.pdf"
        path "sex_assessment_plot.pdf"
        path "metadata_with_sex_info.tsv"
        path "tmm.log"

    script:
        """
        module load miniforge/24.3.0-py3.11
        source activate base
        
        # Run script
        python ${params.scripts_dir}/tmm_norm_pca_sex.py \\
            --metadata ${meta_file} \\
            --gene_counts ${gene_count_file} \\
            --gtf ${gtf_file} \\
            --mappability ${mapp_file} \\
            --output_norm norm_tmm.tsv \\
            --output_pca pca_tmm.tsv \\
            --output_plot_pdf pca_tmm.pdf \\
            --output_sex_plot_pdf sex_assessment_plot.pdf \\
            --xist_threshold 10.0 \\
            --rps4y1_threshold 1.0 \\
            --low_expression_threshold 0.1 \\
            --sample_expression_frac 0.2 \\
            --affected_expression_frac 0.2
            # --exclude_unclear_sex # Include if you want to filter samples with unclear sex.
            # --skip_mappability_filter # Include if you want to implement mappability filtration.
        
        cp .command.log tmm.log
        """
}
