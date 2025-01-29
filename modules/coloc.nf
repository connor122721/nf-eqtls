#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Process for analyzing eqtl output
process prepGWAS {

    // Publish the output to the specified directory
    shell = '/usr/bin/env bash'
    publishDir "${params.out}/gwas", mode: 'copy'
    memory = '30 GB'
    threads = 4

    output:
        path("processed*")

    script:
        """
        module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1

        # Run analysis script
        Rscript ${params.scripts_dir}/prep_GWAS_eQTL_for_coloc.R \\
            --gwas ${params.gwas}
        """
}

// Run Coloc
process runColoc {

    // Publish the output to the specified directory
    shell = '/usr/bin/env bash'
    publishDir "${params.out}/coloc", mode: 'copy'

    input:
        tuple path(eqtl),
        val(chromosome),
        val(best_k)

    output:
        path("*")

    script:
        """
        module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1

        # Get input files
        processed_GWAS=${params.out}/gwas/*${chromosome}.rds
        firstRun=${params.out}/tensorqtl/*${chromosome}_MaxPC30.cis_qtl.txt.gz

        # Run coloc analysis script
        Rscript ${params.scripts_dir}/run_cisNominaleQTL_coloc.R \\
            --gwas \${processed_GWAS} \\
            --eqtl ${eqtl} \\
            --shortList \${firstRun} \\
            --chromosome ${chromosome} \\
            --N_gwas 1665481 \\
            --N_eqtl 516
        """
}
