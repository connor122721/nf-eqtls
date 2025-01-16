#!/usr/bin/env nextflow
/*
 * main_DNA_mega.nf
 * A Nextflow DSL2 pipeline that:
 *  1. Reads a list of chromosomes.
 *  2. Concatenates all non-thinned VCFs.
 *  3. Run King on extracted full VCF.
 *  4. For each chromosome, extracts and indexes a VCF using bcftools.
 *  5. Calculates linkage disequilibrium and prunes the VCF using plink2.
 *  6. Thins VCF using plink2 --bp-space 250.
 *  7. Concatenates all thinned VCFs into one final VCF.
 *  8. Create GDS object from final thinned VCF.
 *  9. Create SNP PCA and output covariates.
 */

nextflow.enable.dsl=2

// ---------------------------
// (A) Define Channels
// ---------------------------
chroms = Channel
    .fromPath(params.chrom_file)
    .splitText()
    .map { it.trim() }
    .filter { it } // remove empty lines

// ---------------------------
// (B) Define Processes
// ---------------------------

/**
 * ExtractAndIndexVCF
 */
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

/**
 * LDPruning
 */
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

        echo "Start LD pruning for chromosome: ${chromosome}"
        date

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

        echo "Finished LD pruning for chromosome: ${chromosome}"
        date
        """
}

/**
 * VCFThin
 * Thins the pruned VCF using plink2 --bp-space 250
 */
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

        echo "Start thinning VCF for chromosome: ${chromosome}"
        date

        plink2 \\
            --threads ${params.threads} \\
            --memory 18000 \\
            --vcf ${pruned_vcf} \\
            --bp-space 250 \\
            --recode vcf bgz \\
            --out ${chromosome}.filt.thin250.TOPchef

        tabix -p vcf ${chromosome}.filt.thin250.TOPchef.vcf.gz

        echo "Finished thinning VCF for chromosome: ${chromosome}"
        date
        """
}

// ---------------------------
// (C) Define Workflow
// ---------------------------

// Side nextflow modules to run
include { ConcatVCF as ConcatVCF_1 } from './modules/concatvcf.nf'
include { ConcatVCF as ConcatVCF_2 } from './modules/concatvcf.nf'
include { King } from './modules/king.nf'
include { VCF_to_GDS } from './modules/SNP_PCA.nf'
include { SNP_PCA } from './modules/SNP_PCA.nf'

// Run Full Workflow! 
workflow {

    // 1) ExtractAndIndexVCF
    def extracted_ch = ExtractAndIndexVCF(chroms)

    // 2) Concat all unthinned VCFs
    def unfilteredList_ch = extracted_ch.map { it[1] }.collect()
    def unfilteredOut_ch  = Channel.value(params.out_vcf_full)
    def unfiltered_Concat = ConcatVCF_1(unfilteredList_ch, unfilteredOut_ch)
    
    // 3) King - kinship
    def kingOut = King(unfiltered_Concat.map { it[0]})

    // 4) LD Pruning
    def pruned_ch = LDPruning(extracted_ch)

    // 5) Thinning
    def thinned_ch = VCFThin(pruned_ch)

    // 6) Concat all thinned VCFs
    def thinnedVcfList_ch = thinned_ch.map { it[1] }.collect()
    def thinnedOutName_ch = Channel.value(params.out_vcf_filt)
    def fin_thinned_vcf = ConcatVCF_2(thinnedVcfList_ch, thinnedOutName_ch)
    
    // 7) Convert VCF to GDS.
    def gds = VCF_to_GDS(fin_thinned_vcf.map { it[0]})

    // 8) SNP-based PCA and plotting.
    def pca = SNP_PCA(gds)

}