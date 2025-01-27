#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Process for analyzing eqtl output
process analysis_eqtl_saturation {

    // Publish the output to the specified directory
    shell = '/usr/bin/env bash'
    publishDir "${params.out}/analyses", mode: 'copy'

    input:
        path(cis_eqtls)

    output:
        path("best_k_eqtls"), emit: bestK
        path("*pdf")
        path("qtl.rna.saturation.rds")

    script:
        """
        module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1

        # Make file for all eQTLs
        ls *cis_qtl.txt.gz > cis_eqtl.list

        # Run analysis script
        Rscript ${params.scripts_dir}/analysis_eQTL_saturation.R \\
            --list_of_eqtls "cis_eqtl.list" \\
            --gtf ${params.gtf_streamlined}
        
        # Finish
        echo "Finish"
        """
}