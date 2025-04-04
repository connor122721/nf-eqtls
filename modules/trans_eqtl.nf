#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// TensorQTL submission process
process TensorTransQTL {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/trans_qtl", mode: 'copy'
    threads = 8
    memory = '30 GB'
    // debug true

    input:
        tuple val(chromosome), 
              val(covariate), 
              val(pc), 
              path(bed_files),
              val(plink_prefix)    

    output:
        tuple path("*cis_qtl.txt.gz"),
            val("${params.out}/tensorqtl"), emit: tensorqtl_out 

    script:
        """
        module load miniforge/24.3.0-py3.11
        source activate qtl

        # Use tensorQTL based on chromosome
        python3 -m tensorqtl \\
            ${params.out}/bedfiles/${plink_prefix} \\
            ${params.out}/eqtl/*.sort.bed \\
            topchef_${chromosome}_MaxPC${pc} \\
            --maf_threshold 0.01 \\
            --covariates ${params.out}/eqtl/${covariate} \\
            --mode trans
        """
}
