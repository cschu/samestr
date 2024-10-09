#!/usr/bin/env python3
import argparse
import os
import pathlib

from collections import Counter
from os.path import basename

import numpy as np

from samestr.convert.buffered_reader import stream_file

BASES = {b: i for b, i in zip('ACGT', range(4))}
STATS_HEADER =  (
	'alignment',
	'genome',
	'mean_cov',
	'median_cov',
	'n_sites',
	'n_gaps',
	'n_covered',
	'n_mono',
	'n_duo',
	'n_tri',
	'n_quat',
	'n_poly',
	'f_covered',
	'f_mono',
	'f_duo',
	'f_tri',
	'f_quat',
	'f_poly',
)


def read_contig_map(contig_map):
    # Contig map
    # cmap[genome] = [contig1, contig2, ...]
    cmap = {}
    # for line in gzip.open(args.map, 'rt'):
    for line in stream_file(contig_map):
        line = line.strip().split()
        genome = line[0]
        contig = line[1]
        cmap.setdefault(genome, []).append(contig)
    return cmap

def initialise_contigs(gene_file):
    # Initialize numpy arrays for each contig
    x = {}
    # for line in gzip.open(args.gene_file, 'rt'):
    for line in stream_file(gene_file):
        line = line.rstrip().split()
        contig = line[0]
        # beg = int(line[2])
        end = int(line[3])
        x[contig] = np.zeros([1, end, 4])
    return x

def add_pileups(positions, contigs):
    # Add kpileup results to numpy arrays
    # with open(args.kp, 'r') as f:
    #     for line in f.readlines():
    # Sample\tContig\tPosition\tGene\tRef\tCodon\tConsensus\tAllele\tCounts\tFrequency
    i = 0
    for sample, contig, pos, gene, ref, codon, consensus, allele, count, freq in positions:
        if allele != "N":
            contigs[contig][i, pos - 1, BASES[allele]] = count

    return contigs

    # for line in stream_file(args.kp):
    #     line = line.rstrip().split()
    #     if len(line) == 10 and line[0] != 'Sample' and line[7] != 'N':
    #         # sample = line[0]
    #         i = 0
    #         contig = line[1]
    #         j = int(line[2])
    #         nt = line[7]
    #         # k = BASES.index(nt)
    #         k = BASES[nt]
    #         count = int(line[8])
    #         x[contig][i, j - 1, k] = count

def read_kpileup(kpileup_file):
    for i, line in enumerate(stream_file(kpileup_file)):
        if i > 0:
            yield line.rstrip().split()
        # line = line.rstrip().split()
        # if len(line) == 10 and line[0] != 'Sample':

def kp2np(kpileups, contig_map, sample, gene_file, output_dir):

    pathlib.Path(output_dir).mkdir(exist_ok=True, parents=True)

    x = initialise_contigs(gene_file)
    x = add_pileups(kpileups, x)
    cmap = read_contig_map(contig_map)

    # with open(os.path.join(output_dir, f"{sample}.aln_stats.txt"), "wt") as stats_out:
        # print(*STATS_HEADER, sep="\t", file=stats_out,)

    # Concatenate contigs
    # -------------------
    m, k = 1, 4
    # y = {
    #     genome: np.zeros([m, sum(np.shape(x[c])[1] for c in contigs), k])
    #     for genome, contigs in cmap.items()
    # }
    y = {}
    # for genome in cmap:
    for genome, contigs in cmap.items():
        # m = 1
        n = sum(np.shape(x[c])[1] for c in contigs)
        # k = 4
        y[genome] = np.zeros([m, n, k])

        # Add alignment data
        beg = 0
        end = 0
        for contig in contigs:
            end += np.shape(x[contig])[1]
            y[genome][0, beg:end, :] = x[contig]
            beg = end

        species = y[genome]
        cov = species.sum(axis=2)

        # coverage [depth]
        # mean_cov = round(np.mean(cov), 4)
        # median_cov = round(np.median(cov), 4)

        # coverage [width]
        n_sites = species.shape[1]
        n_gaps = (cov == 0).sum()
        n_covered = n_sites - n_gaps

        # n of variant sites, monomorphic, .., polymorphic
        # n_mono = ((species > 0).sum(axis=2) == 1).sum()
        # n_duo = ((species > 0).sum(axis=2) == 2).sum()
        # n_tri = ((species > 0).sum(axis=2) == 3).sum()
        # n_quat = ((species > 0).sum(axis=2) == 4).sum()
        # n_poly = ((species > 0).sum(axis=2) > 1).sum()

        # fraction of covered sites,
        # fraction of covered sites with variant, monomorphic, .., polymorphic
        # if not n_covered == 0:
        # if n_covered != 0:
        #     f_covered = round(n_covered / n_sites, 4)
        #     f_mono = round(n_mono / n_covered, 4)
        #     f_duo = round(n_duo / n_covered, 4)
        #     f_tri = round(n_tri / n_covered, 4)
        #     f_quat = round(n_quat / n_covered, 4)
        #     f_poly = round(n_poly / n_covered, 4)

        # else:
        #     f_covered, f_mono, f_duo, \
        #         f_tri, f_quat, f_poly = 0, 0, 0, 0, 0, 0

        # print(
        #     f"{genome}.{sample}",
        #     genome, mean_cov, median_cov, n_sites, n_gaps,
        #     n_covered, n_mono, n_duo, n_tri, n_quat, n_poly, f_covered, f_mono,
        #     f_duo, f_tri, f_quat, f_poly,
        #     sep="\t", file=stats_out,
        # )
        
        # only write to numpy file if there is coverage left after convert criteria
        # if not n_covered == 0:
        # if n_covered != 0:
        if n_covered:
            np.savez_compressed(
                os.path.join(output_dir, f"{genome}.{sample}"),
                y[genome],
                allow_pickle=True
            )




def main():
    # Input arguments
    # ---------------

    parser = argparse.ArgumentParser()
    parser.add_argument('--kp', help='Kpileup alignments (.kp.txt)')
    parser.add_argument('--map', help='Map of genomes to contigs (tab-delimited)')
    parser.add_argument('--sample', help='Sample Name', type=str, required=True)
    parser.add_argument('--gene-file', help='kpileup gene file')
    parser.add_argument('--output-dir', help='Output dir', default='./')
    args = parser.parse_args()

    kp2np(read_kpileup(args.kp), args.map, args.sample, args.gene_file, args.output_dir)

    # Read data
    # ---------

    # Numpy alignments
    # x[contig] = numpy alignment
    # y[genome] = concatenated numpy alignments
    # Merge kpileups from multiple samples. Write dictionary of (M, N, 4) numpy arrays where:
    # M = samples
    # N = alignment sites
    # 4 = nucleotides (ACGT)
    # Entry (i,j,k) of this array corresponds to the count of nucleotide k at position j of sample i

    # Initialize data
    # nts = 'ACGT'
    # M = 1

    # # Initialize numpy arrays for each contig
    # x = {}
    # # for line in gzip.open(args.gene_file, 'rt'):
    # for line in stream_file(args.gene_file):
    #     line = line.rstrip().split()
    #     contig = line[0]
    #     # beg = int(line[2])
    #     end = int(line[3])
    #     x[contig] = np.zeros([1, end, 4])
    
    
    # # Add kpileup results to numpy arrays
    # # with open(args.kp, 'r') as f:
    # #     for line in f.readlines():
    # for line in stream_file(args.kp):
    #     line = line.rstrip().split()
    #     if len(line) == 10 and line[0] != 'Sample' and line[7] != 'N':
    #         sample = line[0]
    #         i = 0
    #         contig = line[1]
    #         j = int(line[2])
    #         nt = line[7]
    #         # k = BASES.index(nt)
    #         k = BASES[nt]
    #         count = int(line[8])
    #         x[contig][i, j - 1, k] = count
    

    # Sample list
    # M = [sample1]
    # M = np.array([args.sample])

    # Contig map
    # cmap[genome] = [contig1, contig2, ...]
    # cmap = {}
    # # for line in gzip.open(args.map, 'rt'):
    # for line in stream_file(args.map):
    #     line = line.strip().split()
    #     genome = line[0]
    #     contig = line[1]
    #     cmap.setdefault(genome, []).append(contig)
    #     # if genome not in cmap:
    #     #     cmap[genome] = []
    #     # cmap[genome].append(contig)

    # Create dir if not exists
    # pathlib.Path(args.output_dir).mkdir(exist_ok=True, parents=True)
    return None
    # Mapping stats
    cols = '\t'.join([
        'alignment', 'genome', 'mean_cov', 'median_cov', 'n_sites', 'n_gaps',
        'n_covered', 'n_mono', 'n_duo', 'n_tri', 'n_quat', 'n_poly', 'f_covered',
        'f_mono', 'f_duo', 'f_tri', 'f_quat', 'f_poly'
    ])
    stats = [cols]

    # Concatenate contigs
    # -------------------
    y = {}
    # for genome in cmap:
    for genome, contigs in cmap.items():
        # contigs = cmap[genome]

        # Initialize array
        # m = len(M)
        m = 1
        n = sum([np.shape(x[c])[1] for c in contigs])
        k = 4
        y[genome] = np.zeros([m, n, k])

        # Add alignment data
        beg = 0
        end = 0
        for contig in contigs:
            end += np.shape(x[contig])[1]
            y[genome][0, beg:end, :] = x[contig]
            beg = end

        species = y[genome]
        cov = species.sum(axis=2)

        # coverage [depth]
        mean_cov = round(np.mean(cov), 4)
        median_cov = round(np.median(cov), 4)

        # coverage [width]
        n_sites = species.shape[1]
        n_gaps = (cov == 0).sum()
        n_covered = n_sites - n_gaps

        # n of variant sites, monomorphic, .., polymorphic
        n_mono = ((species > 0).sum(axis=2) == 1).sum()
        n_duo = ((species > 0).sum(axis=2) == 2).sum()
        n_tri = ((species > 0).sum(axis=2) == 3).sum()
        n_quat = ((species > 0).sum(axis=2) == 4).sum()
        n_poly = ((species > 0).sum(axis=2) > 1).sum()

        # fraction of covered sites,
        # fraction of covered sites with variant, monomorphic, .., polymorphic
        # if not n_covered == 0:
        if n_covered != 0:
            f_covered = round(n_covered / n_sites, 4)
            f_mono = round(n_mono / n_covered, 4)
            f_duo = round(n_duo / n_covered, 4)
            f_tri = round(n_tri / n_covered, 4)
            f_quat = round(n_quat / n_covered, 4)
            f_poly = round(n_poly / n_covered, 4)

        else:
            f_covered, f_mono, f_duo, \
                f_tri, f_quat, f_poly = 0, 0, 0, 0, 0, 0

        np_filepath = '%s/%s.%s' % (args.output_dir, genome, args.sample)
        stat = [
            basename(np_filepath), genome, mean_cov, median_cov, n_sites, n_gaps,
            n_covered, n_mono, n_duo, n_tri, n_quat, n_poly, f_covered, f_mono,
            f_duo, f_tri, f_quat, f_poly
        ]

        # always write to sample stats
        stat = [str(s) for s in stat]
        stats.append('\t'.join(stat))

        # only write to numpy file if there is coverage left after convert criteria
        # if not n_covered == 0:
        if n_covered != 0:
            np.savez_compressed(np_filepath, y[genome], allow_pickle=True)


    with open('%s/%s.aln_stats.txt' % (args.output_dir, args.sample), 'w') as file:
        file.write('\n'.join(stats))


if __name__ == "__main__":
    main()