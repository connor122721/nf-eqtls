#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// --------------------------------------------------
//  Process for reformatting gtf for QTL mapping
// --------------------------------------------------

// Process for reformatting gtf
process reformat_tss_gtf {

    // Publish the output to the specified directory
    shell = '/usr/bin/env bash'
    publishDir "${params.out}/genome_files", mode: 'copy'

    input:
        path raw_gtf

    output:
        path "${params.gtf_name}.bed", emit: gtf

    script:
        """
        module load miniforge/24.3.0-py3.11
        source activate base

        # Collapse annotations to one per transcript
        python ${params.scripts_dir}/collapse_annotations.py \\
            ${raw_gtf} \\
            collapsed.gtf

        # Run python script to reformat the gtf file
        python ${params.scripts_dir}/reformat_TSS_gtf.py \\
            --gtf collapsed.gtf \\
            --output ${params.gtf_name}.bed
        """
}
