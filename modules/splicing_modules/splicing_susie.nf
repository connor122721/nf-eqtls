#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// TensorQTL submission process
process TensorQTLSusie {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/cis_sqtl_susie", mode: 'copy'
    threads = 8
    memory = '30 GB'

    input:
        tuple val(chromosome), 
              val(covariate), 
              val(pc),
              val(plink_prefix),
              path(nominal_files)

    output:
        path("*SuSiE*")

    script:
        """
        module load miniforge/24.3.0-py3.11
        source activate qtl

        # Use tensorQTL based on chromosome
        python3 -m tensorqtl \\
            ${params.out}/bedfiles/${plink_prefix} \\
            ${params.out}/splicing/cluster/*leafcutter.bed.gz \\
            topchefSplice_${chromosome}_MaxPC${pc} \\
            --maf_threshold 0.01 \\
            --covariates ${params.out}/splicing/cluster/${covariate} \\
            --cis_output ${params.out}/splicing/sqtls_nominal/topchefSplice_${chromosome}_MaxPC${pc}.cis_qtl_pairs.${chromosome}.parquet \\
            --mode cis_susie
        """
}