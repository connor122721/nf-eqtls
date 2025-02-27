#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// TensorQTL submission process
process TensorQTLSubmission {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/tensorqtl", mode: 'copy'
    threads = 8
    memory = '15 GB'
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
            --mode cis
        """
}

// TensorQTL submission process for nominal p-value
// These are very large output files so we will only use a subset
process TensorQTLNominal {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/tensorqtl_nominal", mode: 'copy'
    threads = 8
    memory = '20 GB'

    input:
        tuple val(chromosome), 
              val(covariate), 
              val(pc), 
              path(bed_files),
              val(plink_prefix)    

    output:
        tuple path("*parquet"),
            val(chromosome)


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
            --mode cis_nominal
        """
}
