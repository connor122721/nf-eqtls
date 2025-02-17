// TensorQTL submission process
process TensorQTLSusie {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/susie", mode: 'copy'
    threads = 8
    memory = '15 GB'

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
        module load miniforge/24.3.0-py3.11
        source activate qtl

        # Use tensorQTL based on chromosome
        python3 -m tensorqtl \\
            ${params.out}/bedfiles/${plink_prefix} \\
            ${params.out}/eqtl/*.sort.bed \\
            topchef_${chromosome}_MaxPC${pc} \\
            --maf_threshold 0.01 \\
            --covariates ${params.out}/eqtl/${covariate} \\
            --cis_output ${params.out}/tensorqtl_nominal/topchef_${chromosome}_MaxPC${pc}.cis_qtl_pairs.${chromosome}.parquet \\
            --mode cis_susie
        """
}