#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Run linkage analyses
process LD_CandidateGenes {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/linkage", mode: 'copy'

    input:
        path(input_vcf)
        val(chrom)
        val(start)
        val(end)
        val(gene)

    output:
        path("*")
  
    script:
        """
        # Modules to load
        module load bcftools
        module load tabix
        module load plink

        # Samples used
        samples=${params.out}/eqtl/topchef_samples_1_15_25.txt

        # Display window information
        echo "Processing Gene: ${gene} (Chromosome: ${chrom})"

        ###### Process and run pairwise LD ######

        # Subset the VCF file to the 1 Mb window
        temp_vcf=temp_${chrom}_${gene}.vcf.gz
        bcftools view -r ${chrom}:${start}-${end} ${input_vcf} -S \${samples} -Oz -o \${temp_vcf}
        tabix -p vcf \${temp_vcf}

        # Generate the input file in plink format
        plink \\
            --memory 18000 \\
            --vcf \${temp_vcf} \\
            --r2 \\
            --ld-window-r2 0.0 \\
            --geno 0.999 \\
            --maf 0.01 \\
            --double-id \\
            --keep-allele-order \\
            --out topchef_${chrom}_${gene}_${start}_${end} \\
            --allow-extra-chr \\
            --set-missing-var-ids @:#\$1,\$2
        """
}