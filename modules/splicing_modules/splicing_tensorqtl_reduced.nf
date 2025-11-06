#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// 1) Run regtools
process runRegtools {

    container 'docker://griffithlab/regtools:release-1.0.0'
    shell = '/usr/bin/env bash'
    publishDir "${params.out}/splicing/junc", mode: 'copy'
    errorStrategy = 'ignore'
    time = '2h'
    memory = '5 GB'
    cpus = 1

    input:
        tuple val(rnaid),
              val(dnaid),
              path(bam),
              path(bamIndex)

    output:
        path("*junc")
  
    script:
        """
        regtools junctions extract \\
            -a 8 -m 50 -M 500000 \\
            -s XS \\
            ${bam} \\
            -o ${dnaid}.junc
        """
}

// 2) Clustering junction files for QTL mapping using GTEx scripts
process runLeafcutterClusterGTEx {

    container 'docker://francois4/leafcutter:latest'
    shell = '/usr/bin/env bash'
    publishDir "${params.out}/splicing/cluster", mode: 'copy'
    errorStrategy = 'ignore'
    time = '5h'
    memory = '50 GB'
    cpus = 4

    input:
        path(juncFiles)

    output:
        path("topchef*")

    script:
        """
        # Create list of junction files
        ls *.junc > junc.txt

        # Run clustering algorithm
        python3 ${params.scripts_dir}/cluster_prepare_fastqtl.py \\
            junc.txt \\
            ${params.out}/genome_files/gencode.v34.GRCh38.genes.exons.txt.gz \\
            ${params.gtf} \\
            "topchef" \\
            ${params.sample_participant_map} \\
            --num_pcs 100 \\
            --min_clu_reads 10 \\
            --min_clu_ratio 0.001 \\
            --leafcutter_dir ${params.scripts_dir}
        """
}

// 3) Reformatting junction files for QTL mapping
process reformatLeaf {

    container 'library://connmurr243/rcoloc/rtidycoloc'
    shell = '/usr/bin/env bash'
    publishDir "${params.out}/splicing/cluster", mode: 'copy'
    errorStrategy = 'ignore'
    memory = '15 GB'
    cpus = 1

    input:
        path("*")

    output:
        tuple path("*splicepc*"),
              val("${params.out}/splicing/cluster"), emit: reformat_out

    script:
        """
        # Create merged data
        Rscript ${params.scripts_dir}/reformat_sqtl.R \\
            --metadata ${params.metadata} \\
            --pca *leafcutter.PCs.txt \\
            --pca_snp ${params.out}/pca/filt_dna_pc*
        """
}

// 4) TensorQTL submission process - saturation analyses for covariate models
process TensorQTLSubmission {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/splicing/sqtls", mode: 'copy'
    cpus = 8
    memory = '10 GB'

    input:
        tuple val(chromosome), 
              val(covariate), 
              val(pc), 
              val(plink_prefix),
              val(reformat_dir)   

    output:
        tuple path("*qtl.txt.gz"),
            val("${params.out}/splicing/sqtls"), emit: tensorqtl_out 

    script:
        """
        module load miniforge/24.3.0-py3.11
        source activate qtl

        # Use tensorQTL based on chromosome
        python3 -m tensorqtl \\
            ${params.out}/bedfiles/${plink_prefix} \\
            ${params.out}/splicing/cluster/*leafcutter.bed.gz \\
            topchefSplice_${chromosome}_MaxPC${pc} \\
            --maf_threshold 0.01 \\
            --covariates ${params.out}/splicing/cluster/${covariate} \\
            --mode cis
        """
}

// 5) TensorQTL submission process for nominal p-value
process TensorQTLNominal {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/splicing/sqtls_nominal", mode: 'copy'
    cpus = 12
    memory = '25 GB'

    input:
        tuple val(chromosome), 
              val(covariate), 
              val(pc), 
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
            ${params.out}/splicing/cluster/*leafcutter.bed.gz \\
            topchefSplice_${chromosome}_MaxPC${pc} \\
            --maf_threshold 0.01 \\
            --covariates ${params.out}/splicing/cluster/${covariate} \\
            --mode cis_nominal
        """
}

// 6) TensorQTL Submission - trans-sQTLs
process TensorTransQTL {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/splicing/trans_sqtl", mode: 'copy'
    threads = 8
    memory = '50 GB'

    input:
        tuple val(chromosome), 
              val(covariate), 
              val(pc),
              val(plink_prefix),
              path(nominal_files)

    output:
        path("*parquet") 

    script:
        """
        module load miniforge/24.3.0-py3.11
        source activate qtl

        # Use tensorQTL based on chromosome
        python3 -m tensorqtl \\
            ${params.out}/bedfiles/${plink_prefix} \\
            ${params.out}/splicing/cluster/*leafcutter.bed.gz \\
            topchefSplice_${chromosome}_MaxPC${pc} \\
            --maf_threshold 0.01 \\
            --covariates ${params.out}/splicing/cluster/${covariate} \\
            --mode trans \\
            --return_r2
        """
}

// Import modules
include { analysis_sqtl_saturation } from '../analysis_eqtl_saturation.nf'
include { runColocSQTL as runColoc_jurgens } from '../coloc.nf'
include { analysisColocSQTL } from '../coloc.nf'
include { TensorQTLSusie } from './splicing_susie.nf'

workflow {

    // Define constants
    best_k = Channel.of("25")  // Or load from file if needed
    pre_jurgens = Channel.of("jurgens24")
    N_jurgens = Channel.of(955733, 516).toList()

    // Chromosome list
    chroms = Channel
        .from((1..22).collect { it.toString() } + ['X'])
        .map { "chr${it}" }

    // Create channel for nominal TensorQTL results
    TensorQTLNominal_out = Channel
        .fromPath("${params.out}/splicing/sqtls_nominal/*MaxPC25*.parquet")
        .map { path ->
            // Extract chromosome name from filename (e.g., "chr3" from "..._chr3_...")
            def fname = path.getName()
            def chrom = (fname =~ /(chr[0-9XY]+)/)[0][1]
            tuple(path, chrom)
        }

    // Check mapping
    //TensorQTLNominal_out.view()

    // Run coloc analysis
    colocJurgens = runColoc_jurgens(
        TensorQTLNominal_out
            .combine(best_k)
            .combine(pre_jurgens)
            .combine(N_jurgens))

    wd3 = colocJurgens.outDir.unique().collect()
    ColocGenes = analysisColocSQTL(wd3.unique())
}