#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// trans-eQTL submission process
process TensorTransQTL {

    container 'library://connmurr243/wgs/topchef_tensorqtl_python.sif:latest'
    shell = '/usr/bin/env bash'
    publishDir "${params.out}/trans_qtl", mode: 'copy'
    threads = 8
    memory = '30 GB'

    input:
        tuple val(chromosome), 
              val(covariate), 
              val(pc), 
              path(bed_files),
              val(plink_prefix)    

    output:
        path("*parquet") 

    script:
        """
        # module load miniforge/24.3.0-py3.11
        # source activate qtl

        # Use tensorQTL based on chromosome
        python3 -m tensorqtl \\
            ${params.out}/bedfiles/${plink_prefix} \\
            ${params.out}/eqtl/*.sort.bed \\
            topchef_${chromosome}_MaxPC${pc} \\
            --maf_threshold 0.01 \\
            --covariates ${params.out}/eqtl/${covariate} \\
            --mode trans \\
            --return_r2
            
        """
}

// trans-eQTL susie submission process
process TensorTransQTLSusie {

    container 'library://connmurr243/wgs/topchef_tensorqtl_python.sif:latest'
    shell = '/usr/bin/env bash'
    publishDir "${params.out}/trans_susie", mode: 'copy'
    threads = 8
    memory = '30 GB'

    input:
        tuple val(chromosome), 
              val(covariate), 
              val(pc), 
              path(bed_files),
              val(plink_prefix)    

    output:
        path("*SuSiE*") 

    script:
        """
        # module load miniforge/24.3.0-py3.11
        # source activate qtl

        # Use tensorQTL based on chromosome
        python3 -m tensorqtl \\
            ${params.out}/bedfiles/${plink_prefix} \\
            ${params.out}/eqtl/*.sort.bed \\
            topchef_${chromosome}_MaxPC${pc} \\
            --maf_threshold 0.01 \\
            --covariates ${params.out}/eqtl/${covariate} \\
            --mode trans_susie \\
            --return_r2
        """
}
