#!/usr/bin/env python3
import argparse
import gzip
import re
import subprocess
import sys

from collections import defaultdict, Counter
from dataclasses import dataclass

from samestr.convert.buffered_reader import stream_file


CIGAR_RE = re.compile(r'(\d+)([MIDNSHP])')


parser = argparse.ArgumentParser(
    description='This script is for parsing the BAM file and look for reads overlapping with the target genes and report the pileup.')
parser.add_argument('sample_id', help='sample ID')
parser.add_argument('bam_file', help='an aligned bam file')
parser.add_argument(
    'gene_file', help='a tab-delimited file of six columns in this order: contigId, geneId, begin, end, strand, DNA transcript seq. (Note: begin<end)')
parser.add_argument('min_bq', type=int,
                    help='the minimum base quality score of a sequenced base')
parser.add_argument('min_mq', type=int,
                    help='the minimum MQ mapping score of the aligned reads')
parser.add_argument('min_d', type=int, help='the minimum depth, an integer.')
args = parser.parse_args()

usage = f"""{sys.argv[0]} <sampleId> <bam file> <gene file> <minimum BQ> <minimum MQ> <minimum D> <maximum D>

{args.sample_id}: sample ID
{args.bam_file}: an aligned bam file
{args.gene_file}: a tab-delimited file of six columns in this order: contigId, geneId, begin, end, strand, DNA transcript seq. (Note: begin<end)
{args.min_bq}: the minimum base quality score of a sequenced base
{args.min_mq}: the minimum MQ mapping score of the aligned reads
{args.min_d}: the minimum depth, an integer.

This script is for parsing the BAM file and look for reads overlapping with the target genes and report the pileup.
Reads must pass the minimum MQ threshold.
For each position, an allele must have a depth >= <minimum D>.
The frequency of each allele is calculated after the filterings.

Note: the quality score must be Sanger Phred score (ascii 33).

Dependencies:
Use Samtools

CIGAR parser: the embedded cigar parser is not versatile, please review it to make sure that it is handling the cigar code appropriately.


The script generates a tab-delimited table directly to STDOUT.
Each row is a base of the queried genes.
A base could be reported repeatedly up to four times as four rows when polymorphisms are observed: A,T, G or C.
The columns are:
sample ID
contig ID
This Base Position
Gene ID
Ref allele of this Base
Condon Position (if coding gene)
Observed Consensus Allele of this Base (in this BAM)
Observed Allele
Coverage Depth of the Observed Allele
Allele Frequency of the Observed Allele

"""

@dataclass
class Gene:
    gene_id: str = None
    contig_id: str = None
    begin: int = 0
    end: int = 0
    strand: str = None
    seq: str = None

def parse_gene(file):
    """
    Parse the input gene file

    Args:
        gene_file (str): the gene file name

    Returns:
        dict: a dictionary containing the gene name as key and the contig, start, end, strand, and sequence as values
    """
    # data = {}
    # with gzip.open(file, "rt") as f:
    #     for line in f:
    for line in stream_file(file):
        line = line.strip()
        contig_id, gene_id, begin, end, strand, seq = line.split("\t")
        # data[gene_id] = Gene(
        yield Gene(
            gene_id,
            contig_id,
            int(begin),
            int(end),
            strand,
            seq,
        )
        # data[gene_id] = {
        #     'contig': contig_id,
        #     'begin': int(begin),
        #     'end': int(end),
        #     'strand': strand,
        #     'seq': seq
        # }
    # return data


def parse_bases(genes):
    """
    Go through each gene and add the nucleotide positions to a dictionary indexed by contigs

    Args:
        genes (dict): a dictionary containing gene data

    Returns:
        dict: a dictionary indexed by contigs and containing the gene name, reference nucleotide, and codon position as values
    """
    # nuc = defaultdict(dict)
    nuc = {}
    # for g, gene_data in genes.items():
    for gene_data in genes:
        # begin = gene_data['begin']
        # end = gene_data['end']
        # c = gene_data['contig']
        # strand = gene_data['strand']
        # temp = list(gene_data['seq'])

        # for i in range(begin, end + 1):
        #     codon_pos = (i - begin + 1) % 3
        #     if strand == '-' and codon_pos != 2:
        #         codon_pos = 1 if codon_pos == 0 else 0
        #     codon_pos = 3 if codon_pos == 0 else codon_pos
        #     nuc[c][i] = f"{g}\t{temp[i-begin]}\t{codon_pos}"

        for i, base in enumerate(gene_data.seq, start=gene_data.begin): 

            # codon_pos = (i - gene_data.begin + 1) % 3
            codon_pos = i % 3
            if gene_data.strand == '-' and codon_pos != 2:
                codon_pos = 1 if codon_pos == 0 else 0
            codon_pos = 3 if codon_pos == 0 else codon_pos
            # nuc[gene_data.contig_id][i] = f"{gene_data.gene_id}\t{base}\t{codon_pos}"
            nuc.setdefault(gene_data.contig_id, {})[i] = (gene_data.gene_id, base, codon_pos)

    return nuc

def decode_cigar(cigar):
    """
    Decode the cigar string

    Args:
        cigar (str): the cigar string

    Returns:
        str: the decoded cigar string
    """
    cigar_parts = CIGAR_RE.findall(cigar)
    # new_string = ''.join(c * int(n) for n, c in cigar_parts)
    # return new_string
    # return [c * int(n) for n, c in cigar_parts]
    return (cc for op in (c * int(n) for n, c in cigar_parts) for cc in op)


def convert_qual(qual_string):
    """
    Convert the quality string to a list of quality scores

    Args:
        qual_string (str): the quality string

    Returns:
        list: a list of quality scores
    """
    # scores = [ord(q) - 33 for q in qual_string]
    # return scores
    return (ord(q) - 33 for q in qual_string)


def pileup(sample_id, bam_file, gene_file, min_bq, min_mq, min_depth):
    """
    Parse the BAM file and look for reads overlapping with the target genes and report the pileup

    Args:
        sample_id (str): the sample ID
        bam_file (str): the BAM file name
        gene_file (str): the gene file name
        min_bq (int): the minimum base quality score of a sequenced base
        min_mq (int): the minimum MQ mapping score of the aligned reads
        min_depth (int): the minimum depth, an integer.

    Returns:
        None
    """
    genes = parse_gene(gene_file)
    bases = parse_bases(genes)

    # f_table = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    f_table = {}

    with subprocess.Popen(["samtools", "view", bam_file], stdout=subprocess.PIPE, universal_newlines=True) as bam_process:
        for line in bam_process.stdout:
            # qname, flag, rname, begin, mapq, cigar, mrnm, mpos, isize, seq, qual, *info = line.strip().split('\t')
            # if rname == "*" or int(mapq) < min_mq or rname not in bases:
            #     continue
            _, _, rname, begin, mapq, cigar, _, _, _, seq, qual, *_ = line.strip().split("\t")

            if int(mapq) < min_mq:
                continue
            
            contig_bases = bases.get(rname)

            if not contig_bases:
                continue

            begin = int(begin)
            # end = begin + len(seq) - 1
            qual_scores = iter(convert_qual(qual))

            aln_string = decode_cigar(cigar)
            seq_bases = iter(seq)
            # qual_scores = iter(qual_scores)
            # ci = s

            # new = []
            # read_i = 0
            p = begin

            for i, cigar_op in enumerate(aln_string, start=begin):
                if cigar_op == "H":
                    continue
                elif cigar_op == "D":
                    base = "-"
                else:
                    # read_i += 1
                    cur_base, cur_qual = next(seq_bases), next(qual_scores)
                    if cigar_op in ("I", "S"):
                        continue
                    # if qual_scores[read_i] < min_bq:
                    if cur_qual < min_bq:
                        base = "-"
                    else:
                        base = cur_base
                # new.append(base)
                p += 1
                is_position_of_interest = contig_bases.get(p)

                if base != "-" and is_position_of_interest is not None:
                    f_table.setdefault(rname, {}).setdefault(p, Counter())[base] += 1


            # for i, nuc in enumerate(new, start=begin):
            #     if nuc != "-" and bases[rname].get(i):
            #         f_table[rname][i][nuc] += 1
            # for cigar_i in range(len(ci)):
            #     base = "-"
            #     if ci[cigar_i] == "D":
            #         base = "-"
            #     elif ci[cigar_i] == "H":
            #         continue
            #     elif ci[cigar_i] in ["I", "S"]:
            #         read_i += 1
            #         continue
            #     elif qual_scores[read_i] < min_bq:
            #         base = "-"
            #         read_i += 1
            #     else:
            #         base = b[read_i]
            #         read_i += 1

            #     new.append(base)

            # b = new
            # for i in range(begin, begin + len(b)):
            #     nuc = b[i - begin]

            #     if bases[rname].get(i):
            #         if nuc != "-":
            #             f_table[rname][i][nuc] += 1

            # for i, nuc in enumerate(new, start=begin):
            #     if nuc != "-" and bases[rname].get(i):
            #         f_table[rname][i][nuc] += 1

    for c in f_table:
        print(f"{c}===")
    print("Sample\tContig\tPosition\tGene\tRef\tCodon\tConsensus\tAllele\tCounts\tFrequency")

    for c, positions in f_table.items():
        # for pos in sorted(f_table[c]):
        for pos, nucs in sorted(positions.items(), key=lambda x:x[0]):
            total = 0
            major = ""

            # for nuc in sorted(f_table[c][pos]):
            for nuc, counts in sorted(nucs.items(), key=lambda x:x[0]):
                # counts = f_table[c][pos][nuc]
                # if counts < min_depth:
                #     continue
                if counts >= min_depth:
                    total += counts
                    if not major or nucs[major] < counts:
                        major = nuc
                    # elif f_table[c][pos][major] < counts:
                    # elif nucs[major] < counts:
                    #     major = nuc

            # for nuc in sorted(f_table[c][pos]):
            for nuc, counts in sorted(nucs.items(), key=lambda x:x[0]):
                # counts = f_table[c][pos][nuc]
                # if counts < min_depth:
                #     continue
                if counts >= min_depth:
                    percent = 100 * counts / total

                    # print(f"{sample_id}\t{c}\t{pos}\t{bases[c][pos]}\t{major}\t{nuc}\t{counts}\t{percent:.0f}")
                    print(sample_id, c, pos, *bases[c][pos], major, nuc, counts, f"{percent:.0f}", sep="\t")


def main():

    pileup(args.sample_id, 
           args.bam_file, 
           args.gene_file, 
           args.min_bq, 
           args.min_mq, 
           args.min_d)

if __name__ == "__main__":
    main()
