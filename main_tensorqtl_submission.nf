#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * A Nextflow pipeline that:
 *  1. Generates combinations of 1 through 30 and every combination of the 23 chromosomes.
 *  2. Creates bed files for each chromosome.
 *  3. Submits TensorQTL jobs for each chromosome and covariate file.
 */

// ---------------------------
// (A) Define Channels
// ---------------------------

chroms = Channel
    .fromPath(params.chrom_file)
    .splitText()
    .map { it.trim() }

pcs = Channel.from(1..30)

// Combine all chromosome / PC 
chrom_covs = chroms
    .combine(pcs)
    .map { [it[0], "topchef_cov_RNApc1_${it[1]}_1.15.25.txt", it[1]] }

// ---------------------------
// (B) Define Processes
// ---------------------------

// Create bed files for each chromosome
process CreateBedFiles {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/bedfiles", mode: 'copy'
    //debug true

    input:
        val chromosome

    output:
        tuple val(chromosome), path("${chromosome}*.bed")

    script:
        """
        module load bcftools
        module load plink/2.00a20230303

        # Convert to bed file using plink
        plink2 \\
            --memory 18000 \\
            --threads ${params.threads} \\
            --mac 1 \\
            --keep ${params.out}/eqtl/topchef_samples*.txt \\
            --vcf ${params.out}/raw/freeze.10b.${chromosome}.pass_only.phased.TOPchef.vcf.gz \\
            --make-bed \\
            --set-all-var-ids @:# \\
            --rm-dup force-first \\
            --output-chr chrM \\
            --allow-extra-chr \\
            --out ${chromosome}_1.17.25.TOPchef
        """
}

// TensorQTL submission process
process TensorQTLSubmission {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/tensorqtl", mode: 'copy'
    debug true

    input:
        tuple val(chromosome), 
              val(covariate), 
              val(pc), 
              path(bedfile)

    output:
        path "*cis_qtl.txt.gz"

    script:
        """
        module load miniforge/24.3.0-py3.11
        source activate qtl

        # Use tensorQTL based on chromosome
        tensorqtl \\
            ${bedfile} \\
            ${params.out}/eqtl/*.sort.bed \\
            topchef_${chromosome}_MaxPC${pc} \\
            --maf_threshold 0.01 \\
            --covariates ${params.out}/eqtl/${covariate} \\
            --mode cis
        """
}

// ---------------------------
// (C) Define Workflow
// ---------------------------

// Run TensorQTL submission for each chromosome and covariate file
workflow {

    // Run bedfiles
    // chroms.view()
    // chrom_covs.view()
    CreateBedFiles(chroms)

    // Submit TensorQTL jobs
    tensorqtl_input_ch = chrom_covs
        .join(CreateBedFiles.out)
        
    //tensorqtl_input_ch.view()
    // Run tensorQTL by chromosome
    TensorQTLSubmission(tensorqtl_input_ch)
}
