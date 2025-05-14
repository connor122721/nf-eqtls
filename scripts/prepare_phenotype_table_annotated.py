#!/usr/bin/env python
# https://github.com/griffithlab/regtools

import sys
import gzip
import numpy as np
import scipy as sc
import pickle
import re
from optparse import OptionParser

from sklearn.decomposition import PCA
from sklearn import preprocessing
from sklearn import linear_model

from scipy.stats import rankdata
from scipy.stats import norm

def qqnorm(x):
    n = len(x)
    a = 3.0/8.0 if n <= 10 else 0.5
    return norm.ppf((rankdata(x) - a) / (n + 1.0 - 2.0 * a))

def stream_table(f, ss=''):
    fc = '#'
    while fc[0] == "#":
        fc = f.readline().strip()
        head = fc.split(ss)

    for ln in f:
        ln = ln.strip().split(ss)
        attr = {}

        for i in range(len(head)):
            try:
                attr[head[i]] = ln[i]
            except:
                break
        yield attr

def get_chromosomes(ratio_file):
    try:
        open(ratio_file)
    except:
        sys.stderr.write("Can't find %s..exiting\n" % (ratio_file))
        return
    sys.stderr.write("Parsing chromosome names...\n")
    chromosomes = set()
    with gzip.open(ratio_file, 'rt') as f:
        f.readline()
        for line in f:
            chromosomes.add(line.split(":")[0])
    return chromosomes

def get_blacklist_chromosomes(chromosome_blacklist_file):
    if chromosome_blacklist_file:
        with open(chromosome_blacklist_file, 'r') as f:
            return f.read().splitlines()
    else:
        return ["X", "Y"]

def parse_gtf_exons(gtf_file):
    exons_by_chr = {}
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            chrom, _, feature_type, start, end, _, _, _, attributes = fields
            if feature_type != 'exon':
                continue
            match = re.search('gene_name "([^"]+)"', attributes)
            gene_name = match.group(1) if match else "NA"
            start, end = int(start), int(end)
            if chrom not in exons_by_chr:
                exons_by_chr[chrom] = []
            exons_by_chr[chrom].append((start, end, gene_name))
    return exons_by_chr

def find_overlapping_gene(chrom, start, end, exons_by_chr):
    if chrom not in exons_by_chr:
        return "NA"
    for exon_start, exon_end, gene_name in exons_by_chr[chrom]:
        if not (end < exon_start or start > exon_end):
            return gene_name
    return "NA"

def main(ratio_file, chroms, blacklist_chroms, pcs=50, exons_by_chr=None):

    dic_pop, fout = {}, {}
    try:
        open(ratio_file)
    except:
        sys.stderr.write("Can't find %s..exiting\n" % (ratio_file))
        return

    sys.stderr.write("Starting...\n")
    for i in chroms:
        fout[i] = open(ratio_file + ".phen_" + i, 'w')
        fout_ave = open(ratio_file + ".ave", 'w')

    valRows, valRowsnn, geneRows = [], [], []
    header = gzip.open(ratio_file, 'rt').readline().split()[1:]

    for i in fout:
        fout[i].write("\t".join(["#Chr", "start", "end", "ID"] + header + ["Gene"]) + '\n')

    for dic in stream_table(gzip.open(ratio_file, 'rt'), ' '):
        chrom = dic['chrom']
        chr_, s, e, clu = chrom.split(":")
        if chr_ in blacklist_chroms:
            continue
        s, e = int(s), int(e)
        gene_name = find_overlapping_gene(chr_, s, e, exons_by_chr) if exons_by_chr else "NA"

        NA_indices, valRow, aveReads = [], [], []
        tmpvalRow = []

        for i, sample in enumerate(header):
            try:
                count = dic[sample]
                num, denom = count.split('/')
                if float(denom) < 1:
                    tmpvalRow.append("NA")
                    NA_indices.append(i)
                else:
                    count = (float(num) + 0.5) / (float(denom) + 0.5)
                    tmpvalRow.append(count)
                    aveReads.append(count)
            except:
                tmpvalRow.append("NA")
                NA_indices.append(i)

        if tmpvalRow.count("NA") > len(tmpvalRow) * 0.4:
            continue

        ave = np.mean(aveReads)

        for c in tmpvalRow:
            valRow.append(ave if c == "NA" else c)

        if np.std(valRow) < 0.005:
            continue

        if len(valRow) > 0:
            fout[chr_].write("\t".join([chr_, str(s), str(e), chrom] + [str(x) for x in valRow] + [gene_name]) + '\n')
            fout_ave.write(" ".join(["%s" % chrom] + [str(min(aveReads)), str(max(aveReads)), str(np.mean(aveReads))]) + '\n')

            valRowsnn.append(valRow)
            valRow = preprocessing.scale(valRow)

            valRows.append(valRow)
            geneRows.append("\t".join([chr_, str(s), str(e), chrom]))
            if len(geneRows) % 1000 == 0:
                sys.stderr.write("Parsed %s introns...\n" % len(geneRows))

    for i in fout:
        fout[i].close()

    matrix = np.array(valRows)
    for i in range(len(matrix[0, :])):
        matrix[:, i] = qqnorm(matrix[:, i])

    fout = {}
    for i in chroms:
        fn = "%s.qqnorm_%s" % (ratio_file, i)
        print("Outputting: " + fn)
        fout[i] = open(fn, 'w')
        fout[i].write("\t".join(['#Chr', 'start', 'end', 'ID'] + header + ['Gene']) + '\n')

    lst = []
    for i in range(len(matrix)):
        chrom, s, e, id_ = geneRows[i].split("\t")
        s, e = int(s), int(e)
        gene_name = find_overlapping_gene(chrom, s, e, exons_by_chr) if exons_by_chr else "NA"
        lst.append((chrom, int(s), "\t".join([chrom, str(s), str(e), id_] + [str(x) for x in matrix[i]] + [gene_name]) + '\n'))

    lst.sort()
    for ln in lst:
        fout[ln[0]].write(ln[2])

    fout_run = open("%s_prepare.sh" % ratio_file, 'w')
    for i in fout:
        fout[i].close()
        fout_run.write("bgzip -f %s.qqnorm_%s\n" % (ratio_file, i))
        fout_run.write("tabix -p bed %s.qqnorm_%s.gz\n" % (ratio_file, i))
    fout_run.close()

    sys.stdout.write("Use `sh %s_prepare.sh' to create index for fastQTL (requires tabix and bgzip).\n" % ratio_file)

    if pcs > 0:
        pcs = min([len(header), pcs])
        pca = PCA(n_components=pcs)
        pca.fit(matrix)
        pca_fn = ratio_file + ".PCs"
        print("Outputting PCs: " + pca_fn)
        pcafile = open(pca_fn, 'w')
        pcafile.write("\t".join(['id'] + header) + '\n')
        pcacomp = list(pca.components_)

        for i in range(len(pcacomp)):
            pcafile.write("\t".join([str(i + 1)] + [str(x) for x in pcacomp[i]]) + '\n')

        pcafile.close()

if __name__ == "__main__":
    parser = OptionParser(usage="usage: %prog [-p num_PCs] input_perind.counts.gz")
    parser.add_option("-p", "--pcs", dest="npcs", default=50, help="number of PCs output")
    parser.add_option("--ChromosomeBlackList", dest="cbl", default="", help="file of blacklisted chromosomes to exclude from analysis")
    parser.add_option("--gtf", dest="gtf_file", default=None, help="GENCODE GTF file for gene annotation")
    (options, args) = parser.parse_args()

    if len(args) == 0:
        sys.stderr.write("Error: no ratio file provided...\n")
        exit(0)

    exons_by_chr = parse_gtf_exons(options.gtf_file) if options.gtf_file else None

    main(args[0], get_chromosomes(args[0]), get_blacklist_chromosomes(options.cbl), int(options.npcs), exons_by_chr)