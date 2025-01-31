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

// ExtractAndIndexVCF - removes indels of length > 50 and minimum allele count of 1.
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

// Side nextflow modules to run
include { reformat_tss_gtf } from './modules/reformat_TSS_gtf.nf'
include { ConcatVCF as ConcatVCF_1 } from './modules/concatvcf.nf'
include { ConcatVCF as ConcatVCF_2 } from './modules/concatvcf.nf'
include { King } from './modules/king.nf'
include { PlotKinship } from './modules/king.nf'
include { VCF_to_GDS } from './modules/SNP_PCA.nf'
include { SNP_PCA } from './modules/SNP_PCA.nf'
include { reformat_eqtl } from './modules/reformat_eqtl.nf'
include { normalize_and_pca; tmm_pipeline } from './modules/mainRNA_flow.nf'

// Run Full Workflow! 
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
    def split = reformat_eqtl(
        file(params.metadata),
        form_gtf.gtf,
        plotKingOut.related_individuals,
        norm_qc.rna_outliers,
        pca.dna_outliers,
        pca_tmm.norm_gene_count,
        pca_tmm.pca_tmm,
        pca.tensorqtl_pca)

}