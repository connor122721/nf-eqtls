#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Run linkage analyses
process LD_CandidateGenes {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/linkage", mode: 'copy'
    threads = 1
    errorStrategy = 'retry'
    memory = '20 GB'

    input:
        tuple path(input_vcf),
            val(chrom),
            val(start),
            val(end),
            val(gene)

    output:
        path("topchef_${chrom}_${gene}_${start}_${end}*ld.gz")
        path("topchef_${chrom}_${gene}_${start}_${end}.snp_names.txt")
        tuple val(chrom), 
            val(gene), emit: gene_combo
  
    script:
        """
        # Modules to load
        module load bcftools
        module load plink/1.90b7.2

        # Samples used
        samples=${params.out}/eqtl/topchef_samples_1_15_25.txt

        # Display window information
        echo "Processing Gene: ${gene} (Chromosome: ${chrom})"

        ###### Process and run pairwise LD ######

        # Subset the VCF file to the 1 Mb window
        temp_vcf="temp_${chrom}_${gene}.vcf.gz"
        
        bcftools view -r ${chrom}:${start}-${end} \\
        ${params.out}/${input_vcf} \\
        -S \${samples} \\
        -Oz -o \${temp_vcf}

        tabix -p vcf \${temp_vcf}
        output_prefix="topchef_${chrom}_${gene}_${start}_${end}"

        # First, convert the VCF to PLINK binary format (BED/BIM/FAM) using your filters.
        plink \\
            --memory 18000 \\
            --vcf \${temp_vcf} \\
            --geno 0.999 \\
            --maf 0.005 \\
            --double-id \\
            --keep-allele-order \\
            --allow-extra-chr \\
            --list-duplicate-vars suppress-first \\
            --make-bed \\
            --out \${output_prefix}

        # Now, use the binary file (BED) to generate the r-square matrix.
        plink \\
            --memory 18000 \\
            --bfile \${output_prefix} \\
            --r square \\
            --geno 0.999 \\
            --maf 0.005 \\
            --keep-allele-order \\
            --double-id \\
            --list-duplicate-vars suppress-first \\
            --allow-extra-chr \\
            --out \${output_prefix}

        # Compress the LD file
        gzip topchef_${chrom}_${gene}_${start}_${end}.ld
        cut -f4 topchef_${chrom}_${gene}_${start}_${end}.bim > \\
            topchef_${chrom}_${gene}_${start}_${end}.snp_names.txt

        # Clean up
        rm \${temp_vcf}*
        """
}

// Make locuscompare plots
process plotLocusCompare {
    
    // Publish the output to the specified directory
    shell = '/usr/bin/env bash'
    publishDir "${params.out}/coloc_locusPlot", mode: 'copy'
    errorStrategy = 'ignore'
    memory = '20 GB'

    input:
        tuple val(chrom),
            val(gwas),
            val(gene),
            val(eqtl)

    output:
        path("*pdf")

    script:
        """
        module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1

        # Run analysis script
        Rscript ${params.scripts_dir}/analysis_eGene_coloc.R \\
            --gwas ${gwas} \\
            --eqtl ${eqtl} \\
            --chromosome ${chrom} \\
            --gene ${gene}
        """
}

// Run TensorQTL submission for each chromosome and covariate file
workflow {

    // X) Get LD matrix for candidate genes
    def vcf = file("${params.out}/freeze.10b.pass_only.snps_indels50_mac1.phased.TOPchef.vcf.gz")

    // LD List
    Channel
        .fromPath("../data/candgenes.ldlist.txt")
        .splitCsv(skip: 1, strip: true, sep: "\t")
        .map { tuple(vcf, it[2], it[3], it[4], it[1]) }
        .set { ld_list }    

    // Calculate LD for each candidate gene
    def ld = LD_CandidateGenes(ld_list)

    // Collect chromosome and gene output
    ld_chrom_gene = ld.gene_combo
    
    // Define GWAS and eQTL inputs
    gwas_list = Channel.from(["shah20_gwas_HF", "levin22_gwas_HF", "jurgens24_gwas_HF"])
    eqtl_input = Channel.from("MaxPC49")
    chrom_list = Channel.from(1..22).map { "chr${it}" }

    // Combine GWAS, eQTL, chromosome, and gene inputs
    chrom_list
        .combine(gwas_list)
        .combine(ld_chrom_gene, by:0)
        .combine(eqtl_input)
        .set { combined_input }

    // Make locuscompare plots
    plotLocusCompare(combined_input)

}