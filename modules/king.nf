// Run King - kinship analyses
process King {
    shell = '/usr/bin/env bash'
    publishDir "${params.out}/king", mode: 'copy'

    // Collect all .vcf.gz paths from VCFThin into a single list
    input:
        path(input_vcf)

    output:
        path("king*")
  
    script:
        """
        module load plink

        # Convert filtered VCF to bed - also filter MAF > 0.01
        plink2 \\
            --memory 18000 \\
            --threads ${params.threads} \\
            --vcf ${input_vcf} \\
            --set-all-var-ids chr@:# \\
            --rm-dup force-first \\
            --keep-allele-order \\
            --allow-extra-chr \\
            --make-bed \\
            --maf 0.01 \\
            --out ${params.king_bed_out}

        # Run king on bed
        ~/king \\
            -b ${params.king_bed_out}.bed \\
            --degree 5 \\
            --related \\
            --prefix king.${params.king_bed_out}
        """
}