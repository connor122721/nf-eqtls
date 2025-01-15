#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// ------------------------------
//  Processes related to DNA PCA
// ------------------------------

process VCF_to_GDS {

    // Publish the output to the specified directory
    shell = '/usr/bin/env bash'
    publishDir "${params.out}", mode: 'copy', pattern: "*.gds"
    
    // Define the input and output
    input:
        path(input_vcf_ch)

    output:
        path("*.gds")

    // Define a unique output name based on input file
    script:
        """
        module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1

        # Run the R script
        Rscript ${params.scripts_dir}/vcf2gds.R \
            ${input_vcf_ch} \
            ${params.gds_filename} \
            ${params.threads}
        """
}

process SNP_PCA {

    // Publish the outputs to the specified PCA output directory
    shell = '/usr/bin/env bash'
    publishDir "${params.out}/pca", mode: 'copy'

    // Define the input
    input:
        path(gds_file)
        path(related_individuals)

    // Define the output
    output:
        path("${params.pca_rds}")
        path("${params.pca_plot}")
        path("${params.tensorqtl_pca}"), emit: tensorqtl_pca

    // Define the script
    script:
        """
        module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1

        # Execute the R script
        Rscript ${params.scripts_dir}/topchef_dna_pca.R \
            ${gds_file} \
            ${params.metadata} \
            ${related_individuals} \
            ${params.pca_rds} \
            ${params.pca_plot} \
            ${params.tensorqtl_pca}
        """
}
