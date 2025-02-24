#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Run susie_coloc.R on coloc and QTL files
process SUSIE_Coloc {

    label 'process_high'
    shell = '/usr/bin/env bash'
    publishDir "${params.out}/susie", mode: 'copy'
    errorStrategy = 'ignore'
    
    input:
        tuple val(chromosome),
            path(coloc_file),
            val(qtl_dir),
            val(qtl_prefix),
            val(ld_dir),
            val(N1)
    
    output:
        path("*txt")
        //path("*pdf")
    
    script:
        """
        module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1

        Rscript ${params.scripts_dir}/susie_coloc.R \\
            --chromosome ${chromosome} \\
            --coloc ${coloc_file} \\
            --qtl_dir ${qtl_dir} \\
            --qtl_prefix ${qtl_prefix} \\
            --ld_dir ${ld_dir} \\
            --N ${N1}
        """
}

// Run Susie functions
workflow {

    // Chromosomes
    Channel.from( (1..22).collect{ "chr${it}" } )
         .set { chrom_list }
    
    // Define coloc, qtl, and ld_dir inputs (customize paths as needed)
    coloc_file = Channel.fromPath("${params.out}/coloc/coloc_eqtl_candidates_full.txt")
    qtl_input_dir = Channel.from("${params.out}/tensorqtl_nominal")
    qtl_prefix = Channel.from("MaxPC49")
    ld_dir = Channel.from("${params.out}/linkage")
    N = Channel.from(516)

    // Combine inputs
    chrom_list
        .combine(coloc_file)
        .combine(qtl_input_dir)
        .combine(qtl_prefix)
        .combine(ld_dir)
        .combine(N)
        .set { susie_inputs }

    // Run Susie on candidate coloc regions
    SUSIE_Coloc(susie_inputs)

}
