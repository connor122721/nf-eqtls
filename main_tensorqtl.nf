#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// ---------------------------
// (A) Define Channels
// ---------------------------

// Chromosomes
chroms = Channel
    .fromPath(params.chrom_file)
    .splitText()
    .map { it.trim() }

// Number of RNA PCs to test
pcs = Channel.from(1..50)

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
            --keep ${params.out}/eqtl/topchef_samples*.txt \\
            --vcf ${params.out}/raw/freeze.10b.${chromosome}.pass_only.phased.TOPchef.vcf.gz \\
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
        tuple path("*parquet"),
            val(chromosome)


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
            --mode cis_nominal
        """
}

// ---------------------------
// (C) Define Workflow
// ---------------------------

include { analysis_eqtl_saturation } from './modules/analysis_eqtl_saturation.nf'
include { prepGWAS as prepGWAS_levin } from './modules/coloc.nf'
include { prepGWAS as prepGWAS_shah } from './modules/coloc.nf'
include { prepGWAS as prepGWAS_jurgens } from './modules/coloc.nf'
include { runColoc as runColoc_levin } from './modules/coloc.nf'
include { runColoc as runColoc_shah } from './modules/coloc.nf'
include { runColoc as runColoc_jurgens } from './modules/coloc.nf'
include { analysisColoc } from './modules/coloc.nf'
//include { LD_CandidateGenes } from './modules/ld.nf'

// Run TensorQTL submission for each chromosome and covariate file
workflow {

    // 1) Make bedfiles
    bed = CreateBedFiles(chroms)
    plink_prefix_ch = bed.plink_prefix

    // Submit TensorQTL jobs
    tensorqtl_input_ch = chrom_covs
        .combine(plink_prefix_ch, by: 0)

    // 2) Run tensorQTL by chromosome
    TensorQTLSubmission(tensorqtl_input_ch)

    // Extract unique path from TensorQTLSubmission
    Channel
        TensorQTLSubmission.out
        .map { tuple -> tuple[0] }
        .collect()
        .set { outi }
    
    // 3) Extract Best-K to run nominal p-value
    analysis_eqtl_saturation(outi)    

    // Extract Best K
    analysis_eqtl_saturation.out.bestK
        .splitText()
        .set { best_k }
    
    // Run tensorQTL - nominal p-value
    chrom_covs
        .combine(plink_prefix_ch, by: 0)
        .filter { it[2] == 49 }  // Filter for PC equal to best k
        .set { tensorqtl_input_nom_ch }

    // 4) Run tensorQTL - nominal p-value
    TensorQTLNominal(tensorqtl_input_nom_ch)

    // 5) Prep GWAS for coloc
    prepGWAS_out_levin = prepGWAS_levin(params.gwas_levin,
                                    "levin22_gwas_HF",
                                    "TRUE")

    prepGWAS_out_shah = prepGWAS_shah(params.gwas_shah,
                                    "shah20_gwas_HF",
                                    "TRUE")

    prepGWAS_out_jurgens = prepGWAS_jurgens(params.gwas_jurgens,
                                    "jurgens24_gwas_HF",
                                    "FALSE")

    // 6) Run coloc analyses
    N_levin = Channel.of(1665481, 516).toList()
    N_shah = Channel.of(977323, 516).toList()
    N_jurgens = Channel.of(955733, 516).toList()
    
    // Run 6a
    colocLevin = runColoc_levin(TensorQTLNominal.out
                   .combine(best_k)
                   .combine(prepGWAS_out_levin.prefix)
                   .combine(N_levin))

    // Run 6b
    colocShah = runColoc_shah(TensorQTLNominal.out
                  .combine(best_k)
                  .combine(prepGWAS_out_shah.prefix)
                  .combine(N_shah))

    // Run 6c
    colocJurgens = runColoc_jurgens(TensorQTLNominal.out
                  .combine(best_k)
                  .combine(prepGWAS_out_jurgens.prefix)
                  .combine(N_jurgens))

    // 7) Get candidate eGenes that are colocalized and prep for LD analysis
    wd1 = colocLevin.outDir.unique().collect()
    wd2 = colocShah.outDir.unique().collect()
    wd3 = colocJurgens.outDir.unique().collect()
 
    ColocGenes = analysisColoc(wd1.join(wd2).join(wd3).unique())

    // X) Get LD matrix for candidate genes
    def vcf = file("params.out/freeze.10b.pass_only.snps_indels50_mac1.phased.TOPchef.vcf.gz")
    genes = ColocGenes.candidate_genes.splitText().toList()

    //genes.view()

    //inputLD = vcf.merge(genes)

    //inputLD.view()

    //LD_CandidateGenes()

}