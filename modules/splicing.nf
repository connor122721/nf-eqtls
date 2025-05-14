#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Run regtools
process runRegtools {

    container 'docker://griffithlab/regtools:release-1.0.0'
    shell = '/usr/bin/env bash'
    publishDir "${params.out}/splicing/junc", mode: 'symlink'
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

// Intron clustering
process runLeafcutterCluster {

    container 'docker://francois4/leafcutter:latest'
    shell = '/usr/bin/env bash'
    publishDir "${params.out}/splicing/cluster", mode: 'symlink'
    memory = '5 GB'
    threads = 1

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
        python ${params.scripts_dir}/leafcutter_cluster_regtools.py \\
            -j ${chrom}.junc.txt \\
            -m 30 \\
            -o ${chrom} \\
            -l 500000
        """
}

// Reformatting junction files for QTL mapping
process reformatLeaf {

    container 'docker://francois4/leafcutter:latest'
    shell = '/usr/bin/env bash'
    publishDir "${params.out}/splicing/sqtl", mode: 'symlink'
    memory = '5 GB'
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
            -p 70
        """
}

// Run splicing analyses for each chromosome
workflow {

    // STEP 1: Splicing quantfication and set up
    
    // Input parameters
    meta_ch = Channel
        .fromPath('../metadata/metadata_10_17_2024_CSM.txt')
        .splitCsv(strip: true, sep: '\t', header: true)
        .map { row -> tuple(row.SAMPLE_ID_NWD.trim(), 
                            row.SAMPLE_ID_TOR.trim())}

    // Sample list
    samples = Channel
        .fromPath('output/eqtl/topchef_samples_1_15_25.txt')
        .splitCsv(strip: true, sep: '\t')
        .map { it[0].toString() } 
        .join(meta_ch, by:0)
        .filter { nwd, tor -> tor in ['TOR855588','TOR795199','TOR487285'] }  

    // Chromosome list
    chroms = Channel
        .from( (1..22).collect { it.toString() } + ['X'] )
        .map { idx -> "chr${idx}" }

    // Cross‑product
    sample_chrom = samples.combine(chroms)

    // Write out bam list
    sample_chrom_bam = sample_chrom.map { rnaID, dnaID, chr ->
        bamPath = "/standard/vol185/TOPMed/TOPCHef/bams/${dnaID}.Aligned.sortedByCoord.out.md.bam"
        baiPath = "/standard/vol185/TOPMed/TOPCHef/bams/${dnaID}.Aligned.sortedByCoord.out.md.bam.bai"
        tuple(rnaID, dnaID, chr, bamPath, baiPath)}

    // Run regtools
    runRegtools(sample_chrom_bam)

    // Collect junction files
    runRegtools.out
      .groupTuple(by: 1)
      .set { reg_grouped }

    // Run Leafcutter clustering
    runLeafcutterCluster(reg_grouped)

    // Reformat junction files for QTL mapping
    reformatLeaf(runLeafcutterCluster.out.collect())

    // STEPs 2: sQTL Mapping 

}