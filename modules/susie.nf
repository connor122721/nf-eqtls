#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process SUSIE_Coloc {
    publishDir "${params.out}/susie", mode: 'copy'
    container "${params.r_container}"  // if using a container with R
    
    input:
        val(chromosome)
        path(coloc_file)
        path(qtl_file)
        path(ld_dir)
    
    output:
        path "susie_${chromosome}.txt"
    
    script:
        """
        module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1
        Rscript ${params.scripts_dir}/susie_coloc.R \\
            --chromosome ${chromosome} \\
            --coloc ${coloc_file} \\
            --qtl ${qtl_file} \\
            --ld_dir ${ld_dir} \\
            --N ${params.susie_N} \\
            --out susie_${chromosome}.txt
        """
}

workflow {
    // Example: Assuming channels for input files are provided in your pipeline
    Channel.from( (1..22).collect{ "chr${it}" } )
         .set { chrom_list }
    
    // Define coloc, qtl, and ld_dir inputs (customize paths as needed)
    coloc_file = file("${params.input_dir}/coloc/coloc_eqtl_levin22_HF_seVar.txt")
    qtl_file   = file("${params.input_dir}/qtl/run_11_12_24/sample.parquet")
    ld_dir     = file("${params.input_dir}/linkage/out")
    
    chrom_list.map { chr -> tuple(chr, coloc_file, qtl_file, ld_dir) }
              .set { susie_inputs }
    
    SUSIE_Coloc(susie_inputs)
}
