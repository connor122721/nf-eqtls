/**
 * ConcatVCF
 * Concatenates all per-chromosome thinned VCFs
 */
process ConcatVCF {

    container 'library://connmurr243/wgs/wgs_common.sif:latest'
    shell = '/usr/bin/env bash'
    publishDir "${params.out}", mode: 'copy'

    // Collect all .vcf.gz paths from VCFThin into a single list
    input:
        path(input_list_vcfs)
        val(vcf_output)

    output:
        path("${vcf_output}*")
  
    script:
        """
        echo "Concatenating the following files:"
        for vcf in ${input_list_vcfs}; do
          echo " \$vcf"
        done

        # Concatenate all files
        bcftools concat \\
            --threads ${params.threads} \\
            -Oz \\
            -o ${vcf_output} \\
            ${input_list_vcfs}

        tabix -p vcf ${vcf_output}
        """
}