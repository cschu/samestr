#!/usr/bin/env python3
import logging
import pathlib
import re
import sqlite3
import sys

import numpy as np

from samestr.convert.buffered_reader import stream_file


LOG = logging.getLogger(__name__)

CIGAR_RE = re.compile(r'(\d+)([MIDNSHP])')

BASES = {"A": 0, "C": 1, "G": 2, "T": 3}


def initialise_contigs(gene_file):
    # Initialize numpy arrays for each contig
    contigs = {}
    for line in stream_file(gene_file):
        contig, _, _, end, _, _ = line.rstrip().split("\t")
        contigs[contig] = np.zeros([1, int(end), 4])
    return contigs

def initialise_contigs_db(clades, db):

    query_placeholders = ",".join("?" * len(clades))

    contigs = {}

    with sqlite3.connect(db) as conn:
        cursor = conn.cursor()
        cursor.execute(f"SELECT clade.name,marker.name,marker.length FROM marker JOIN clade ON marker.clade_id = clade.id WHERE clade.name IN ({query_placeholders}) ORDER BY clade.name,contig.name", clades)
        for clade, contig, length in cursor.fetchall():
            contigs[contig] = clade, np.zeros([1, length, 4])
    
    return contigs





def decode_cigar(cigar):
    for n, c in CIGAR_RE.findall(cigar):
        if c != "H":
            yield int(n), c

def pileup(bam_stream, gene_file, min_bq, min_mq, min_depth, clades, db, outstream=sys.stdout):
    # contigs = initialise_contigs(gene_file)
    contigs = initialise_contigs_db(clades, db)

    cur_rname = None
    contig = None

    for line in bam_stream:
        _, _, rname, begin, mapq, cigar, _, _, _, seq, qual, *_ = line.strip().split("\t")

        if int(mapq) < min_mq:
            continue

        if rname != cur_rname:
            cur_rname = rname
            _, contig = contigs.get(rname)

        if contig is None:
            continue

        p = int(begin)
        
        bases_and_quals = zip(seq, qual)
        
        for oplen, cigar_op in decode_cigar(cigar):
            if cigar_op == "H":
                continue
            elif cigar_op != "D":
                is_insertion = cigar_op in ("I", "S")
                for pp in range(p, p + oplen):
                    cur_base, cur_qual = next(bases_and_quals)
                    if not is_insertion:
        
                        cur_base = BASES.get(cur_base)

                        if cur_base is not None and (ord(cur_qual) - 33) >= min_bq and pp <= contig.shape[1]:
                            contig[0, pp - 1, cur_base] += 1

                if is_insertion:
                    continue
            
            p += oplen

    for contig_name, (clade, contig) in contigs.items():
        contig[contig < min_depth] = 0
        yield clade, contig_name, contig
    # for c, positions in f_table.items():
    #     for pos, nucs in sorted(positions.items(), key=lambda x:x[0]):
    #         for nuc, counts in sorted(nucs.items(), key=lambda x:x[0]):
    #             if counts >= min_depth:
    #                 print(c, pos, nuc, counts, sep="\t", file=outstream)
    #                 yield c, pos, nuc, counts


def read_contig_map2(contig_map):
    # Contig map
    # cmap[genome] = [contig1, contig2, ...]
    cmap = {}
    for line in stream_file(contig_map):
        line = line.strip().split()
        genome = line[0]
        contig = line[1]
        cmap.setdefault(genome, []).append(contig)
    return cmap

def read_contig_map(contig_map):
    cur_genome, contigs = None, []
    for line in stream_file(contig_map):
        genome, contig = line.rstrip().split("\t")
        if cur_genome != genome:
            if cur_genome is not None and contigs:
                yield cur_genome, contigs
                contigs.clear()
            cur_genome = genome
        contigs.append(contig)
    if cur_genome is not None and contigs:
        yield cur_genome, contigs
        

def process_genome(sample, genome, contig_pileups, output_dir):
    n = sum(np.shape(c[c])[1] for c in contig_pileups)
    y = np.zeros([1, n, 4])

    # Add alignment data
    beg, end = 0, 0
    for contig in contig_pileups:
        end += np.shape(contig)[1]
        y[0, beg:end, :] = contig
        beg = end


    # only write to numpy file if there is coverage left after convert criteria
    cov = y.sum(axis=2)
    n_gaps = (cov == 0).sum()
    n_covered = n - n_gaps

    if n_covered:
        np.savez_compressed(output_dir / f"{genome}.{sample}", y, allow_pickle=True)


def kp2np(kpileups, contig_map, sample, gene_file, output_dir):

    output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    # x = {contig_id: pileups for contig_id, pileups in kpileups}
    # cmap = read_contig_map(contig_map)

    # Concatenate contigs
    # -------------------
    # for genome, contigs in cmap.items():
    # for genome, contigs in read_contig_map(contig_map):

    genome_contigs = []
    cur_genome = None
    for genome, contig, pileups in kpileups:
        if genome != cur_genome:
            if genome_contigs:
                process_genome(sample, cur_genome, genome_contigs, output_dir)
                genome_contigs.clear()
            cur_genome = genome
        genome_contigs.append(contig)
    
    if genome_contigs:
        process_genome(sample, cur_genome, genome_contigs, output_dir)



        # n = sum(np.shape(x[c])[1] for c in contigs)
        # y = np.zeros([1, n, 4])

        # # Add alignment data
        # beg, end = 0, 0
        # for contig in contigs:
        #     end += np.shape(x[contig])[1]
        #     y[0, beg:end, :] = x[contig]
        #     beg = end


        # # only write to numpy file if there is coverage left after convert criteria
        # cov = y.sum(axis=2)
        # n_gaps = (cov == 0).sum()
        # n_covered = n - n_gaps

        # if n_covered:
        #     np.savez_compressed(output_dir / f"{genome}.{sample}", y, allow_pickle=True)
