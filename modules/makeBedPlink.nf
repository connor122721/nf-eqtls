#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Create bed files for each chromosome
process CreateBedFiles {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/bedfiles", mode: 'copy'
    //debug true

    input:
        val chromosome
        path samples

    output:
        tuple val(chromosome),
            path("${chromosome}*"), 
            val("${chromosome}_1.17.25.TOPchef"), emit: plink_prefix

    script:
        """
        module load bcftools
        module load plink/2.00a20230303

        # Convert to bed file using plink
        plink2 \\
            --memory 18000 \\
            --threads ${params.threads} \\
            --mac 1 \\
            --keep ${samples} \\
            --vcf ${params.out}/raw/*.${chromosome}.*vcf.gz \\
            --make-bed \\
            --set-all-var-ids @:# \\
            --rm-dup force-first \\
            --output-chr chrM \\
            --allow-extra-chr \\
            --out ${chromosome}_1.17.25.TOPchef
        
        # Get rid of intermediary log
        rm *TOPchef.log
        """
}
