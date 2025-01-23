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
        tuple val(chromosome),
            path("${chromosome}*.bed"), 
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
    // debug true

    input:
        tuple val(chromosome), 
              val(covariate), 
              val(pc), 
              path(bed_files),
              val(plink_prefix)    

    output:
        tuple path("*cis_qtl.txt.gz"),
            val("${params.out}/tensorqtl"), emit: tensorqtl_out 

    script:
        """
        module load miniforge/24.3.0-py3.11
        source activate qtl

        # Use tensorQTL based on chromosome
        python3 -m tensorqtl \\
            ${params.out}/bedfiles/${plink_prefix} \\
            ${params.out}/eqtl/*.sort.bed \\
            topchef_${chromosome}_MaxPC${pc} \\
            --maf_threshold 0.01 \\
            --covariates ${params.out}/eqtl/${covariate} \\
            --mode cis
        """
}

// TensorQTL submission process for nominal p-value
// These are very large output files so we will only use a subset
process TensorQTLNominal {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/tensorqtl_nominal", mode: 'copy'
    // debug true

    input:
        tuple val(chromosome), 
              val(covariate), 
              val(pc), 
              path(bed_files),
              val(plink_prefix)    

    output:
        path "*cis_qtl.txt.gz"

    script:
        """
        module load miniforge/24.3.0-py3.11
        source activate qtl

        # Use tensorQTL based on chromosome
        python3 -m tensorqtl \\
            ${params.out}/bedfiles/${plink_prefix} \\
            ${params.out}/eqtl/*.sort.bed \\
            topchef_${chromosome}_MaxPC${pc} \\
            --maf_threshold 0.01 \\
            --covariates ${params.out}/eqtl/${covariate} \\
            --mode cis-nominal
        """
}

// ---------------------------
// (C) Define Workflow
// ---------------------------

include { analysis_eqtl_saturation } from './modules/analysis_eqtl_saturation.nf'

// Run TensorQTL submission for each chromosome and covariate file
workflow {

    // Run bedfiles
    // chroms.view()
    // chrom_covs.view()
    bed = CreateBedFiles(chroms)
    plink_prefix_ch = bed.plink_prefix

    // Debugging outputs
    // plink_prefix_ch.view()

    // Submit TensorQTL jobs
    tensorqtl_input_ch = chrom_covs
        .combine(plink_prefix_ch, by: 0)

    // chrom_covs.combine(plink_prefix_ch, by:0).view()
    // Run tensorQTL by chromosome
    TensorQTLSubmission(tensorqtl_input_ch)

    // Run tensorQTL - nominal p-value
    tensorqtl_input_nom_ch = chrom_covs
        .combine(plink_prefix_ch, by: 0)
        .filter { it[2] == 15 }  // Add this line to filter for PC equal to 15

    // tensorqtl_input_nom_ch.view()

    // Extract unique path from TensorQTLSubmission
    Channel
        TensorQTLSubmission.out
        .map { tuple -> tuple[1] }
        .unique()
        .collect()
        .map { list -> list[0] }
        .set { out }

    // analysis_eqtl_saturation(out)

    // Extract Best K to run nominal p-value
    

    // Run tensorQTL - nominal p-value
    // TensorQTLNominal(tensorqtl_input_nom_ch)
}