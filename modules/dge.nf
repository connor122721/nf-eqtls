// Run Differential expression analysis
process dge_deseq {

    label 'process_high'
    shell = '/usr/bin/env bash'
    publishDir "${params.out}/dge", mode: 'copy'

    input:
        path metadata
        path gene_counts
        path gtf
        path outdir

    output:
        path 'sig_dge_DCM.csv'

    script:
        """
        module load miniforge/24.3.0-py3.11
        source activate base

        python ${params.scripts_dir}/dge_deseq2.py \
            --metadata ${metadata} \\
            --gene_counts ${gene_counts} \\
            --gtf ${gtf} \\
            --outdir ${outdir}
        """
}