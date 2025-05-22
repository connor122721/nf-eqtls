#!/usr/bin/env python
# Edited by: Connor Murray (5/14/2025)
# Updated: 5/15/2025 to assign genes to clusters and extract TSS with corrected TSS start/end and phenotype_id containing original intron info

import sys
import gzip
import numpy as np
import re
from optparse import OptionParser
from collections import defaultdict

from sklearn.decomposition import PCA
from sklearn import preprocessing
from scipy.stats import rankdata, norm

from qtl.io import gtf_to_tss_bed

def qqnorm(x):
    n = len(x)
    a = 3.0 / 8.0 if n <= 10 else 0.5
    return norm.ppf((rankdata(x) - a) / (n + 1.0 - 2.0 * a))

def stream_table(f, ss=''):
    fc = '#'
    while fc[0] == "#":
        fc = f.readline().strip()
        head = fc.split(ss)
    for ln in f:
        ln = ln.strip().split(ss)
        attr = {head[i]: ln[i] for i in range(min(len(ln), len(head)))}
        yield attr

def get_chromosomes(ratio_file):
    with gzip.open(ratio_file, 'rt') as f:
        f.readline()
        return set(line.split(":")[0] for line in f)

def get_blacklist_chromosomes(chromosome_blacklist_file):
    if chromosome_blacklist_file:
        with open(chromosome_blacklist_file, 'r') as f:
            return f.read().splitlines()
    return ["X", "Y"]

def parse_gene_intervals(gtf_file):
    gene_intervals = defaultdict(list)
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            chrom, _, feature_type, start, end, _, _, _, attributes = fields
            if feature_type != 'gene':
                continue
            match = re.search(r'gene_id "([^"]+)"', attributes)
            if match:
                gene_id = match.group(1)
                gene_intervals[chrom].append((int(start), int(end), gene_id))
    return gene_intervals

def find_overlapping_gene(chrom, start, end, gene_intervals):
    for gene_start, gene_end, gene_id in gene_intervals.get(chrom, []):
        if not (end < gene_start or start > gene_end):
            return gene_id
    return "NA"

def load_tss_from_gtf(gtf_file):
    tss_df = gtf_to_tss_bed(gtf_file)
    tss_dict = defaultdict(dict)
    for _, row in tss_df.iterrows():
        chrom, start, end, gene_id = row["chr"], int(row["start"]), int(row["end"]), row["gene_id"]
        gene_id = row.name
        if gene_id != "NA":
            tss_dict[chrom][gene_id] = (start, end, gene_id)
    return tss_dict

def main(ratio_file, chroms, blacklist_chroms, pcs=50, tss_dict=None, gene_intervals=None):
    dic_pop, fout = {}, {}
    header = gzip.open(ratio_file, 'rt').readline().split()[1:]

    for chrom in chroms:
        fout[chrom] = open(f"{ratio_file}.phen_{chrom}", 'w')
        fout[chrom].write("\t".join(["#Chr", "start", "end", "phenotype_id"] + header) + '\n')

    fout_ave = open(ratio_file + ".ave", 'w')
    valRows, valRowsnn, geneRows = [], [], []

    for dic in stream_table(gzip.open(ratio_file, 'rt'), ' '):
        chrom = dic['chrom']
        chr_, s, e, cluster_id = chrom.split(":")
        if chr_ in blacklist_chroms:
            continue

        s, e = int(s), int(e)
        gene_id = find_overlapping_gene(chr_, s, e, gene_intervals)
        if gene_id == "NA" or tss_dict is None:
            continue

        tss_info = tss_dict.get(chr_, {}).get(gene_id)
        if not tss_info:
            continue

        tss_start, tss_end, gene_id = tss_info
        phenotype_id = f"{chrom}:{gene_id}"  # cluster string plus gene_id

        values, NA_indices, aveReads = [], [], []
        for i, sample in enumerate(header):
            try:
                count = dic[sample]
                num, denom = count.split('/')
                denom = float(denom)
                if denom < 1:
                    values.append("NA")
                    NA_indices.append(i)
                else:
                    value = (float(num) + 0.5) / (denom + 0.5)
                    values.append(value)
                    aveReads.append(value)
            except:
                values.append("NA")
                NA_indices.append(i)

        if values.count("NA") > len(values) * 0.4:
            continue

        ave = np.mean(aveReads)
        valRow = [ave if v == "NA" else v for v in values]
        if np.std(valRow) < 0.005:
            continue

        fout[chr_].write("\t".join([chr_, str(tss_start), str(tss_end), phenotype_id] + [str(x) for x in valRow]) + '\n')
        fout_ave.write(" ".join([cluster_id] + [str(min(aveReads)), str(max(aveReads)), str(np.mean(aveReads))]) + '\n')

        valRowsnn.append(valRow)
        valRow = preprocessing.scale(valRow)
        valRows.append(valRow)
        geneRows.append("\t".join([chr_, str(tss_start), str(tss_end), phenotype_id]))
        if len(geneRows) % 1000 == 0:
            sys.stderr.write(f"Parsed {len(geneRows)} introns...\n")

    for f in fout.values():
        f.close()

    if len(valRows) == 0:
        print("No valid introns after filtering. Exiting.")
        return

    matrix = np.array(valRows)
    for i in range(len(matrix[0])):
        matrix[:, i] = qqnorm(matrix[:, i])

    fout = {}
    for chrom in chroms:
        fn = f"{ratio_file}.qqnorm_{chrom}"
        fout[chrom] = open(fn, 'w')
        fout[chrom].write("\t".join(['#Chr', 'start', 'end', 'phenotype_id'] + header) + '\n')

    lst = []
    for i in range(len(matrix)):
        chrom, s, e, phenotype_id = geneRows[i].split("\t")
        lst.append((chrom, int(s), "\t".join([chrom, s, e, phenotype_id] + [str(x) for x in matrix[i]]) + '\n'))

    lst.sort()
    for ln in lst:
        fout[ln[0]].write(ln[2])
        fout[ln[0]].flush()

    with open(f"{ratio_file}_prepare.sh", 'w') as fout_run:
        for chrom in fout:
            fout[chrom].close()
            fout_run.write(f"bgzip -f {ratio_file}.qqnorm_{chrom}\n")
            fout_run.write(f"tabix -p bed {ratio_file}.qqnorm_{chrom}.gz\n")

    if pcs > 0:
        pcs = min(len(header), pcs)
        pca = PCA(n_components=pcs)
        pca.fit(matrix)
        with open(ratio_file + ".PCs", 'w') as pcafile:
            pcafile.write("\t".join(['id'] + header) + '\n')
            for i, component in enumerate(pca.components_):
                pcafile.write("\t".join([str(i + 1)] + [str(x) for x in component]) + '\n')

if __name__ == "__main__":
    parser = OptionParser(usage="usage: %prog [-p num_PCs] input_perind.counts.gz")
    parser.add_option("-p", "--pcs", dest="npcs", default=50, help="number of PCs output")
    parser.add_option("--ChromosomeBlackList", dest="cbl", default="", help="file of blacklisted chromosomes")
    parser.add_option("--gtf", dest="gtf_file", default=None, help="GENCODE GTF file for gene annotation")
    (options, args) = parser.parse_args()

    if len(args) == 0:
        sys.stderr.write("Error: no ratio file provided...\n")
        exit(1)

    gene_intervals = parse_gene_intervals(options.gtf_file)
    tss_dict = load_tss_from_gtf(options.gtf_file)
    main(args[0], get_chromosomes(args[0]), get_blacklist_chromosomes(options.cbl), int(options.npcs), tss_dict, gene_intervals)