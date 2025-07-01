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
            --num_pcs 10 \\
            --min_clu_reads 5 \\
            --min_clu_ratio 0.000001 \\
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
    memory = '30 GB'

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
// These are very large output files so we will only use a subset of the saturation tests
process TensorQTLNominal {

    shell = '/usr/bin/env bash'
    publishDir "${params.out}/splicing/sqtls_nominal", mode: 'copy'
    cpus = 8
    memory = '30 GB'

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

// Other scripts for alternative analyses
// Intron clustering
process runLeafcutterCluster {

    container 'docker://francois4/leafcutter:latest'
    shell = '/usr/bin/env bash'
    publishDir "${params.out}/splicing/cluster", mode: 'copy'
    errorStrategy = 'ignore'
    time = '10h'
    memory = '15 GB'
    threads = 4

    input:
        tuple path(juncFiles),
              val(chrom)

    output:
        path("*perind.counts.gz")

    script:
        """
        # Create list of junction files
        ls *${chrom}.junc > ${chrom}.junc.txt

        # Run clustering algorithm
        python3 ${params.scripts_dir}/leafcutter_cluster_regtools.py \\
            -j ${chrom}.junc.txt \\
            -m 30 \\
            -o ${chrom} \\
            -l 500000
        """
}

// Reformatting junction files for QTL mapping
process OLD_reformatLeaf {

    container 'docker://francois4/leafcutter:latest'
    shell = '/usr/bin/env bash'
    publishDir "${params.out}/splicing/sqtl", mode: 'copy'
    errorStrategy = 'ignore'
    memory = '10 GB'
    threads = 1

    input:
        path(juncFiles)

    output:
        path("*qqnorm_chr*")
        path("*counts.gz.PCs")
        path("*phen_chr*")

    script:
        """
        # Create merged junction files
        zcat chr1_perind.counts.gz | head -n 1 > all_clusters.counts
        sed -i 's/.chr1//g' all_clusters.counts # Replace Chromosome from sample names
        for file in chr*_perind.counts.gz; do
            zcat "\${file}" | tail -n +2 >> all_clusters.counts
        done
        gzip all_clusters.counts

        # Run clustering algorithm
        python3 ${params.scripts_dir}/prepare_phenotype_table_annotated.py \\
            all_clusters.counts.gz \\
            --gtf "/standard/vol185/cphg_Manichaikul/users/csm6hg/genome_files/gencode.v34.GRCh38.genes.collapsed_only.gtf" \\
            -p 100
        """
}

// Regtools split by chromosome
process runRegtoolsChrom {

    container 'docker://griffithlab/regtools:release-1.0.0'
    shell = '/usr/bin/env bash'
    publishDir "${params.out}/splicing/junc", mode: 'copy'
    errorStrategy = 'ignore'
    time = '3h'
    memory = '2 GB'
    threads = 1

    input:
        tuple val(rnaid),
              val(dnaid),
              val(chrom),
              path(bam),
              path(bamIndex)

    output:
        tuple path("*junc"),
              val(chrom)
  
    script:
        """
        regtools junctions extract \\
            -a 8 -m 50 -M 500000 \\
            -s XS \\
            -r ${chrom} \\
            ${bam} \\
            -o ${dnaid}.${chrom}.junc
        """
}

// Import modules
include { analysis_sqtl_saturation } from '../analysis_eqtl_saturation.nf'
include { runColocSQTL as runColoc_levin } from '../coloc.nf'
include { runColocSQTL as runColoc_shah } from '../coloc.nf'
include { runColocSQTL as runColoc_jurgens } from '../coloc.nf'
include { analysisColocSQTL } from '../coloc.nf'

// Run splicing analyses for each chromosome
workflow {

    // STEP I: Splicing quantfication and set up
    
    // Input parameters
    meta_ch = Channel
        .fromPath('../metadata/metadata_10_17_2024_CSM.txt')
        //.fromPath('../metadata/metadata_6_17_2025_CSM_reduced.txt')
        .splitCsv(strip: true, sep: '\t', header: true)
        .map { row -> tuple(row.SAMPLE_ID_NWD.trim(), 
                            row.SAMPLE_ID_TOR.trim())}

    // Sample list
    samples = Channel
        .fromPath('output/eqtl/topchef_samples_1_15_25.txt')
        .splitCsv(strip: true, sep: '\t')
        .map { it[0].toString() } 
        .join(meta_ch, by:0)

    // Chromosome list
    chroms = Channel
        .from( (1..22).collect { it.toString() } + ['X'] )
        .map { idx -> "chr${idx}" }

    // Write out bam list
    sample_bam = samples.map { rnaID, dnaID ->
        bamPath = "/standard/vol185/TOPMed/TOPCHef/bams/${dnaID}.Aligned.sortedByCoord.out.md.bam"
        baiPath = "/standard/vol185/TOPMed/TOPCHef/bams/${dnaID}.Aligned.sortedByCoord.out.md.bam.bai"
        tuple(rnaID, dnaID, bamPath, baiPath)}

    // Run regtools
    runRegtools(sample_bam)

    // Collect junction files
    Channel 
        runRegtools.out
        .collect()
        .set { reg_grouped }

    // Run Leafcutter clustering
    runLeafcutterClusterGTEx(reg_grouped)

    // Reformat junction files for QTL mapping
    reformatLeaf(runLeafcutterClusterGTEx.out)

    // STEPs II: sQTL Mapping 

    // Number of Splicing PCs to test
    pcs = Channel.from(1..100)

    // Combine all chromosome / PC 
    chrom_covs = chroms
        .combine(pcs)
        .combine(reformatLeaf.out.reformat_out.collect())
        .map { [it[0], "topchef_cov_splicepc1_${it[1]}_06.23.25.txt", it[1], "${it[0]}_1.17.25.TOPchef", it[2]] }

    // Submit TensorQTL - saturation QTL jobs
    TensorQTLSubmission(chrom_covs)

    // Extract unique path from TensorQTLSubmission
    Channel
        TensorQTLSubmission.out
        .map { tuple -> tuple[0] }
        .collect()
        .set { outi }
    
    // Extract Best-K to run nominal p-value
    analysis_sqtl_saturation(outi)    

     // Extract Best K
    analysis_sqtl_saturation.out.bestK
        .splitText()
        .set { best_k }

    // Submit TensorQTL nominal jobs - for best model
    //chrom_covs
    //    .filter { it[2] == 2 }  // Filter for PC equal to best k
    //    .set { tensorqtl_input_nom_ch }

    // Run tensorQTL - cis-nominal p-value
    //TensorQTLNominal(tensorqtl_input_nom_ch)

    // Step III: Colocalization 
    
    // Run coloc analyses
    //N_levin = Channel.of(1665481, 516).toList()
    //N_shah = Channel.of(977323, 516).toList()
    //N_jurgens = Channel.of(955733, 516).toList()
    //pre_lev = Channel.of("levin22")
    //pre_shah = Channel.of("shah20")
    //pre_jurgens = Channel.of("jurgens24")

    // Run a
    //colocLevin = runColoc_levin(TensorQTLNominal.out
    //               .combine(best_k)
    //              .combine(pre_lev)
    //               .combine(N_levin))

    // Run b
    //colocShah = runColoc_shah(TensorQTLNominal.out
    //              .combine(best_k)
    //              .combine(pre_shah)
    //             .combine(N_shah))

    // Run c
    //colocJurgens = runColoc_jurgens(TensorQTLNominal.out
    //              .combine(best_k)
    //              .combine(pre_jurgens)
    //              .combine(N_jurgens))

    // Get candidate eGenes that are colocalized and prep for LD analysis
    //wd1 = colocLevin.outDir.unique().collect()
    //wd2 = colocShah.outDir.unique().collect()
    //wd3 = colocJurgens.outDir.unique().collect() 
    //ColocGenes = analysisColocSQTL(wd1.join(wd2).join(wd3).unique())
}