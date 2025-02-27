#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * A Nextflow pipeline that:
 *  Reformats GTF for transcription start sites.
 *  Extracts and indexes a VCF using bcftools.
 *  Concatenates all non-thinned VCFs.
 *  Runs King on extracted full VCF.
 *  Calculates linkage disequilibrium and prunes the VCF.
 *  Thins VCF.
 *  Concatenates all thinned VCFs into one final VCF.
 *  Creates GDS object from final thinned VCF.
 *  Creates SNP PCA and outputs covariates.
 *  Reformats PC covariates and VCF for cis-eQTL pipeline.
 *  Normalizes RNA counts.
 *  Performs PCA on normalized RNA counts.
 *  Runs eQTL pipeline.
 *  Runs colocalization analysis.
 */

// ---------------------------
// (A) Define Channels
// ---------------------------
chroms = Channel
    .fromPath(params.chrom_file)
    .splitText()
    .map { it.trim() }
    .filter { it }

// ---------------------------
// (B) Define Processes
// ---------------------------

// Removes indels of length > 50 and minimum allele count of 1.
process ExtractAndIndexVCF {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/raw", mode: 'copy'

    input:
        val(chromosome)

    output:
        tuple val(chromosome), 
              path("freeze.10b.${chromosome}.pass_only.phased.TOPchef.vcf.gz"),
              path("freeze.10b.${chromosome}.pass_only.phased.TOPchef.vcf.gz.tbi")

    script:
        """
        module load bcftools
        module load htslib

        echo "Generating new file for chromosome: ${chromosome}"

        bcftools view \\
            "${params.wd}/freeze.10b.${chromosome}.pass_only.phased.bcf" \\
            --types snps,indels \\
            -i 'TYPE="snp" || (TYPE="indel" && ILEN<50)' \\
            --min-ac 1 \\
            -r ${chromosome} \\
            --threads "${params.threads}" \\
            --samples-file "${params.samp}" \\
            -Oz -o freeze.10b.${chromosome}.pass_only.phased.TOPchef.vcf.gz

        tabix -p vcf freeze.10b.${chromosome}.pass_only.phased.TOPchef.vcf.gz

        echo "Done with chromosome: ${chromosome}"
        """
}

// LDPruning
process LDPruning {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/ldPrune", mode: 'copy'

    input:
        tuple val(chromosome), 
              path(input_vcf), 
              path(input_vcf_index)

    output:
        tuple val(chromosome),
              path("freeze.10b.${chromosome}.filt.TOPchef.vcf.gz"),
              path("freeze.10b.${chromosome}.filt.TOPchef.vcf.gz.tbi")

    script:
        """
        module load plink/2.00a20230303
        module load htslib

        # First pass: generate prune.in / prune.out
        plink2 \\
            --memory 18000 \\
            --threads ${params.threads} \\
            --vcf ${input_vcf} \\
            --set-all-var-ids chr@:# \\
            --rm-dup force-first \\
            --keep-allele-order \\
            --indep-pairwise 100 10 0.1 \\
            --allow-extra-chr \\
            --out freeze.10b.${chromosome}.pass_only.snps_indels50_mac1_phased.TOPchef

        # Second pass: remove pruned sites (keep only prune.in)
        plink2 \\
            --memory 18000 \\
            --threads ${params.threads} \\
            --vcf ${input_vcf} \\
            --set-all-var-ids chr@:# \\
            --rm-dup force-first \\
            --keep-allele-order \\
            --allow-extra-chr \\
            --extract freeze.10b.${chromosome}.pass_only.snps_indels50_mac1_phased.TOPchef.prune.in \\
            --recode vcf bgz \\
            --out freeze.10b.${chromosome}.filt.TOPchef

        tabix -p vcf freeze.10b.${chromosome}.filt.TOPchef.vcf.gz
        """
}

// Thins the pruned VCF using plink2 --bp-space 250
process VCFThin {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/thin", mode: 'copy'

    input:
        tuple val(chromosome),
              path(pruned_vcf),
              path(pruned_vcf_index)

    output:
        tuple val(chromosome),
              path("${chromosome}.filt.thin250.TOPchef.vcf.gz"),
              path("${chromosome}.filt.thin250.TOPchef.vcf.gz.tbi")

    script:
        """
        module load plink/2.00a20230303
        module load htslib

        plink2 \\
            --memory 18000 \\
            --threads ${params.threads} \\
            --vcf ${pruned_vcf} \\
            --bp-space 250 \\
            --recode vcf bgz \\
            --out ${chromosome}.filt.thin250.TOPchef

        tabix -p vcf ${chromosome}.filt.thin250.TOPchef.vcf.gz
        """
}

// ---------------------------
// (C) Define Workflow
// ---------------------------

// Modules to run for QC and processing
include { reformat_tss_gtf } from './modules/reformat_TSS_gtf.nf'
include { ConcatVCF as ConcatVCF_1 } from './modules/concatvcf.nf'
include { ConcatVCF as ConcatVCF_2 } from './modules/concatvcf.nf'
include { King } from './modules/king.nf'
include { PlotKinship } from './modules/king.nf'
include { VCF_to_GDS } from './modules/SNP_PCA.nf'
include { SNP_PCA } from './modules/SNP_PCA.nf'
include { reformat_eqtl } from './modules/reformat_eqtl.nf'
include { normalize_and_pca; tmm_pipeline } from './modules/mainRNA_flow.nf'

// Modules to run for eQTL mapping and colocalization
include { CreateBedFiles } from './modules/makeBedPlink.nf'
include { TensorQTLSubmission } from './modules/tensorqtl.nf'
include { TensorQTLNominal } from './modules/tensorqtl.nf'

include { analysis_eqtl_saturation } from './modules/analysis_eqtl_saturation.nf'
include { prepGWAS as prepGWAS_levin } from './modules/coloc.nf'
include { prepGWAS as prepGWAS_shah } from './modules/coloc.nf'
include { prepGWAS as prepGWAS_jurgens } from './modules/coloc.nf'

include { coloc as runColoc_levin } from './modules/coloc.nf'
include { coloc as runColoc_shah } from './modules/coloc.nf'
include { coloc as runColoc_jurgens } from './modules/coloc.nf'
include { analysisColoc } from './modules/coloc.nf'

// Run Full Workflow! Woo! 
workflow {

    // 1) ExtractAndIndexVCF
    def extracted_ch = ExtractAndIndexVCF(chroms)
    
    // 2) Reformatted GTF
    def form_gtf = reformat_tss_gtf(params.gtf)
    
    // 3) Concat all unthinned VCFs
    def unfilteredList_ch = extracted_ch.map { it[1] }.collect()
    def unfilteredOut_ch  = Channel.value(params.out_vcf_full)
    def unfiltered_Concat = ConcatVCF_1(unfilteredList_ch, unfilteredOut_ch)
    
    // 4) King - kinship
    def kingOut = King(unfiltered_Concat.map { it[0]})

    // 5) Plot kinship and identify related individuals
    def plotKingOut = PlotKinship(kingOut)

    // 6) LD Pruning
    def pruned_ch = LDPruning(extracted_ch)

    // 7) Thinning
    def thinned_ch = VCFThin(pruned_ch)

    // 8) Concat all thinned VCFs
    def thinnedVcfList_ch = thinned_ch.map { it[1] }.collect()
    def thinnedOutName_ch = Channel.value(params.out_vcf_filt)
    def fin_thinned_vcf = ConcatVCF_2(thinnedVcfList_ch, thinnedOutName_ch)
    
    // 9) Convert VCF to GDS.
    def gds = VCF_to_GDS(fin_thinned_vcf.map { it[0]})

    // 10) SNP-based PCA and plotting.
    def pca = SNP_PCA(gds, plotKingOut.related_individuals)

    // 11) Normalize RNA counts - drops samples.
    def norm_qc = normalize_and_pca(
        file(params.metadata), 
        file(params.mapp_file),
        form_gtf.gtf,
        file(params.gene_count_file))

    // 12) PCA on normalized RNA counts - drops genes.
    def pca_tmm = tmm_pipeline(
        file(params.metadata), 
        file(params.mapp_file), 
        form_gtf.gtf, 
        file(params.gene_count_file))

    // 13) Reformat PC covariates and VCF for cis-eQTL pipeline.
    def reform = reformat_eqtl(
        file(params.metadata),
        form_gtf.gtf,
        plotKingOut.related_individuals,
        norm_qc.rna_outliers,
        pca.dna_outliers,
        pca_tmm.norm_gene_count,
        pca_tmm.pca_tmm,
        pca.tensorqtl_pca)

    // Move onto eQTL pipeline below!

    // 1) Make bedfiles
    bed = CreateBedFiles(chroms, reform.sample_list)
    plink_prefix_ch = bed.plink_prefix

    // Number of RNA PCs to test
    pcs = Channel.from(1..50)

    // Combine all chromosome / PC 
    chrom_covs = chroms
        .combine(pcs)
        .map { [it[0], "topchef_cov_RNApc1_${it[1]}_1.15.25.txt", it[1]] }

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

}