#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Process for analyzing eqtl output
process prepGWAS {

    // Publish the output to the specified directory
    shell = '/usr/bin/env bash'
    publishDir "${params.out}/gwas", mode: 'copy'
    memory = '30 GB' // This is a memory intensive job!
    threads = 4

    input:
        path(gwas)
        val(prefix)
        val(liftover)

    output:
        path("processed*")
        val(prefix), emit: prefix

    script:
        """
        module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1

        # Run analysis script
        Rscript ${params.scripts_dir}/prep_GWAS_eQTL_for_coloc.R \\
            --gwas ${gwas} \\
            --prefix ${prefix} \\
            --liftover ${liftover}
        """
}

// Run Coloc
process runColoc {

    // Publish the output to the specified directory
    shell = '/usr/bin/env bash'
    publishDir "${params.out}/coloc", mode: 'copy'
    errorStrategy = 'ignore'

    input:
        tuple path(eqtl),
        val(chromosome),
        val(best_k),
        val(gwas_pre),
        val(N_gwas),
        val(N_eqtl)

    output:
        path("*")
        val("${params.out}/coloc"), emit: outDir

    script:
        """
        module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1

        # Get input files
        processed_GWAS=${params.out}/gwas/*${gwas_pre}*${chromosome}.rds
        firstRun=${params.out}/tensorqtl/*${chromosome}_MaxPC49.cis_qtl.txt.gz

        # Run coloc analysis script - Levin 2022 GWAS
        Rscript ${params.scripts_dir}/run_cisNominaleQTL_coloc.R \\
            --gwas \${processed_GWAS} \\
            --eqtl ${eqtl} \\
            --shortList \${firstRun} \\
            --chromosome ${chromosome} \\
            --N_gwas ${N_gwas} \\
            --N_eqtl ${N_eqtl} \\
            --prefix ${gwas_pre}

        echo "Finish!"
        """
}

// Process for analyzing eqtl output
process analysisColoc {

    // Publish the output to the specified directory
    shell = '/usr/bin/env bash'
    publishDir "${params.out}/coloc", mode: 'copy'

    input:
        path(coloc_results)

    output:
        path("coloc_eqtl_candidates_regions.txt"), emit: candidate_genes
        path("*pdf")
        path("coloc_eqtl_candidates_full.txt")

    script:
        """
        module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1

        # Run analysis script
        Rscript ${params.scripts_dir}/analysis_coloc.R \\
            --wd ${coloc_results}
        """
}
