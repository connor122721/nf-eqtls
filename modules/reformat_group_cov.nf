#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// --------------------------------------------------
//  Processes related to covariate reformatting
// --------------------------------------------------

// Process for reformatting the covariates
process reformat_group_covariates {

    // Publish the output to the specified directory
    shell = '/usr/bin/env bash'
    publishDir "${params.out}/bed_group", mode: 'copy'
    
    // Define the input and output
    input:
        tuple val(chrom), 
        val(group),
        path(input_vcf)

    output:
        path("*bed")
        path("*bim")
        path("*.fam")

    script:
        """
        module load plink/2.00a20230303
        
        echo "Start reformatting VCF for chromosome: ${chrom}"
        echo "Group: ${group}"

        # Conditional statement to determine group
        if [ "${group}" == "affected" ]; then
            samp_file="${params.samp_aff}"
        else
            samp_file="${params.samp_non}"
        fi

        echo "Using sample file: \${samp_file}"

        plink2 \\
            --memory 18000 \\
            --threads ${params.threads} \\
            --vcf ${input_vcf} \\
            --make-bed \\
            --set-all-var-ids @:# \\
            --rm-dup force-first \\
            --output-chr chrM \\
            --allow-extra-chr \\
            --keep \${samp_file} \\
            --out ${chrom}.${params.group_bed_out}.${group}
        """
}

// Process for reformatting the covariates
process reformat_covariates {

    // Publish the output to the specified directory
    shell = '/usr/bin/env bash'
    publishDir "${params.out}/bed", mode: 'copy'
    
    // Define the input and output
    input:
        tuple val(chrom),
        path(input_vcf)

    output:
        path("*bed")
        path("*bim")
        path("*.fam")

    script:
        """
        module load plink/2.00a20230303
        
        echo "Start reformatting VCF for chromosome: ${chrom}"

        plink2 \\
            --memory 18000 \\
            --threads ${params.threads} \\
            --vcf ${input_vcf} \\
            --make-bed \\
            --set-all-var-ids @:# \\
            --rm-dup force-first \\
            --output-chr chrM \\
            --allow-extra-chr \\
            --out ${chrom}.${params.group_bed_out}
        """
}