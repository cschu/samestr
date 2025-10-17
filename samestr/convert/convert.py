#!/usr/bin/env python3
import logging
import pathlib
import re
import sys

from collections import defaultdict

import numpy as np

from samestr.convert.buffered_reader import stream_file


LOG = logging.getLogger(__name__)

CIGAR_RE = re.compile(r'(\d+)([MIDNSHP])')

BASES = {"A": 0, "C": 1, "G": 2, "T": 3}


def parse_positions(f):
    positions = {}
    for line in stream_file(f):
        contig_id, _, start, end, _, _ = line.rstrip().split("\t")
        positions[contig_id] = int(start), int(end)
    return positions

def decode_cigar(cigar):
    for n, c in CIGAR_RE.findall(cigar):
        if c != "H":
            yield int(n), c

def convert_qual(qual_string):
    for q in qual_string:
        yield ord(q) - 33


def pileup(bam_stream, gene_file, min_bq, min_mq, min_depth, outstream=sys.stdout):
    # bases = parse_positions(gene_file)
    contigs = initialise_contigs(gene_file)

    # f_table = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    # refs = {}
    # freqs = []
    cur_rname = None

    contig = None
    for line in bam_stream:
        _, _, rname, begin, mapq, cigar, _, _, _, seq, qual, *_ = line.strip().split("\t")

        if int(mapq) < min_mq:
            continue

        if rname != cur_rname:
            cur_rname = rname
            contig = contigs.get(rname)
            
            # contig_bases = bases.get(rname)
            # if contig_bases is not None:
            #     start, end = contig_bases
            #     rindex = refs.setdefault(rname, len(refs))
            #     if rindex >= len(freqs):
            #         freqs.append(np.zeros([1, end, 4]))
            #     # f_table.setdefault(rname, np.zeros([1, end, 4]))

        # if contig_bases is None:
        if contig is None:
            continue
        # start, end = contig_bases

        p = int(begin)
        
        # bases_and_quals = iter(
        #     zip(
        #         (BASES.get(b) for b in seq),
        #         convert_qual(qual)
        #     )
        # )
        bases_and_quals = zip(seq, qual)
        
        for oplen, cigar_op in decode_cigar(cigar):
            if cigar_op == "H":
                continue
            # elif cigar_op == "D":
            #     base = "-"                
            # else:
            elif cigar_op != "D":
                is_insertion = cigar_op in ("I", "S")
                for pp in range(p, p + oplen):
                    cur_base, cur_qual = next(bases_and_quals)
                    if not is_insertion:
                        # base = (cur_base, "-")[(cur_qual < min_bq)]
                        # if base != "-" and start <= pp <= end:
                        # if cur_qual >= min_bq and cur_base is not None and 1 <= pp <= contig.shape[1]:

                        cur_base = BASES.get(cur_base)

                        if cur_base is not None and (ord(cur_qual) - 33) >= min_bq and pp <= contig.shape[1]:

                            # f_table[rname][pp][base] += 1
                            contig[0, pp - 1, cur_base] += 1
                            # f_table[rname][0, pp - 1, cur_base] += 1
                
                if is_insertion:
                    continue
            
            p += oplen

    for contig_name, contig in contigs.items():
        contig[contig < min_depth] = 0
        yield contig_name, contig
    # for c, positions in f_table.items():
    #     for pos, nucs in sorted(positions.items(), key=lambda x:x[0]):
    #         for nuc, counts in sorted(nucs.items(), key=lambda x:x[0]):
    #             if counts >= min_depth:
    #                 print(c, pos, nuc, counts, sep="\t", file=outstream)
    #                 yield c, pos, nuc, counts



def add_pileups(positions, contigs):
    for contig, pos, allele, count in positions:
        if allele != "N":
            contigs[contig][0, pos - 1, BASES[allele]] = count


def read_contig_map(contig_map):
    # Contig map
    # cmap[genome] = [contig1, contig2, ...]
    cmap = {}
    for line in stream_file(contig_map):
        line = line.strip().split()
        genome = line[0]
        contig = line[1]
        cmap.setdefault(genome, []).append(contig)
    return cmap

def initialise_contigs(gene_file):
    # Initialize numpy arrays for each contig
    contigs = {}
    for line in stream_file(gene_file):
        line = line.rstrip().split()
        contig = line[0]
        end = int(line[3])
        contigs[contig] = np.zeros([1, end, 4])
    return contigs



def kp2np(kpileups, contig_map, sample, gene_file, output_dir):

    output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    # x = initialise_contigs(gene_file)
    # add_pileups(kpileups, x)
    x = {contig_id: pileups for contig_id, pileups in kpileups}
    cmap = read_contig_map(contig_map)

    # Concatenate contigs
    # -------------------
    m, k = 1, 4
    # y = {}
    for genome, contigs in cmap.items():
        n = sum(np.shape(x[c])[1] for c in contigs)
        # y[genome] = np.zeros([m, n, k])
        y = np.zeros([m, n, k])

        # Add alignment data
        beg = 0
        end = 0
        for contig in contigs:
            end += np.shape(x[contig])[1]
            # y[genome][0, beg:end, :] = x[contig]
            y[0, beg:end, :] = x[contig]
            beg = end

        # species = y[genome]
        # cov = species.sum(axis=2)
        cov = y.sum(axis=2)

        # n_sites = species.shape[1]
        n_gaps = (cov == 0).sum()
        # n_covered = n_sites - n_gaps
        n_covered = n - n_gaps

        # only write to numpy file if there is coverage left after convert criteria
        if n_covered:
            np.savez_compressed(
                output_dir / f"{genome}.{sample}",
                # y[genome],
                y,
                allow_pickle=True
            )
