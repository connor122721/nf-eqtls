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
        file "norm_medrat.tsv"
        file "pca_medrat.tsv"
        file "pca_medrat_plot.pdf"

    script:
        """
        module load miniforge/24.3.0-py3.11
        source activate base
        
        python ${params.scripts_dir}/medratio_norm_pca.py \\
            --metadata $meta_file \\
            --mappability $mapp_file \\
            --gtf $gtf_file \\
            --gene_counts $gene_count_file \\
            --output_normalized norm_medrat.tsv \\
            --output_pca pca_medrat.tsv \\
            --output_outliers ${params.dna_outliers}
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
        file "norm_tmm.tsv", emit: norm_gene_count
        file "pca_tmm.tsv"
        file "pca_tmm.pdf"
        file "sex_assessment_plot.pdf"

    script:
        """
        module load miniforge/24.3.0-py3.11
        source activate base

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
          # --exclude_unclear_sex   # Uncomment if you want to exclude samples with unclear sex
        """
}
