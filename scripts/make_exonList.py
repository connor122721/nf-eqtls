import pandas as pd
import qtl.annotation

annot = qtl.annotation.Annotation('/standard/vol185/cphg_Manichaikul/users/csm6hg/genome_files/gencode.v34.GRCh38.annotation.gtf')
exon_df = pd.DataFrame([[g.chr, e.start_pos, e.end_pos, g.strand, g.id, g.name]
                        for g in annot.genes for e in g.transcripts[0].exons],
                       columns=['chr', 'start', 'end', 'strand', 'gene_id', 'gene_name'])
exon_df.to_csv('/standard/vol185/cphg_Manichaikul/users/csm6hg/nextflow_dna/output/genome_files/gencode.v34.GRCh38.genes.exons.txt.gz', sep='\t', index=False)