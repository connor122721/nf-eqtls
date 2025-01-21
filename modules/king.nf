#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// ---------------------------------------
//  Processes related to kinship analyses
// ----------------------------------------

// Run King - kinship analyses
process King {
    shell = '/usr/bin/env bash'
    publishDir "${params.out}/king", mode: 'copy'

    // Collect all .vcf.gz paths into a single list
    input:
        path(input_vcf)

    output:
        path("*kin"), emit: king_out
  
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

// Run R script for plotting kinship and identifying related individuals
process PlotKinship {
    shell = '/usr/bin/env bash'
    publishDir "${params.out}/king", mode: 'copy'

    input:
        path(king_files)

    output:
        path("kinship.pdf")
        path("relatedIndividuals.txt"), emit: related_individuals

    script:
        """
        module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1

        Rscript ${params.scripts_dir}/plot_king.R \\
            --king_input ${king_files} \\
            --metadata_file ${params.metadata} \\
            --output_plot kinship.pdf \\
            --output_related relatedIndividuals.txt
        """
}