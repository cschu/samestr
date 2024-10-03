#!/usr/bin/env python3
import re
import sys

from contextlib import nullcontext

from samestr.convert.buffered_reader import stream_file

""" Filter SAM file by percent identity and length """


COMPLEMENT_TRANSLATION = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
CIGAR_RE = re.compile(r'([a-zA-Z])')
HARDCLIP_RE = re.compile(r'H')

def decode_cigar(cigar):
    hbeg = 0
    alnlen = 0
    hend = 0
    tlen = 0
    tabs = CIGAR_RE.split(cigar)[:-1]
    it_tabs = iter(tabs)
    for i, c in enumerate(it_tabs):
        ci = int(c)
        li = next(it_tabs)
        if li == 'H':
            if i == 0:
                hbeg += ci
            else:
                hend += ci
        if li == 'M':
            alnlen += ci
        if li != 'H' and li != 'D':
            tlen += ci

    return hbeg, alnlen, hend, tlen


def reverse_complement(seq):
    return seq.translate(COMPLEMENT_TRANSLATION)[::-1]

def main():
    # initialize variables
    cquery = ''
    cseq = ''
    cqual = ''
    cstrand = ''    

    # read command line arguments and/or stdin
    if (len(sys.argv[1:]) == 2):
        # fhandle = sys.stdin.readlines()
        instream = sys.stdin
        instream_context = nullcontext()
        pctid = float(sys.argv[1])
        minlen = int(sys.argv[2])
    else:
        # instream = instream_context = open(sys.argv[1], 'rt')
        instream = stream_file(sys.argv[1])
        instream_context = nullcontext()
        # fhandle = f.readlines()
        pctid = float(sys.argv[2])
        minlen = int(sys.argv[3])

    id_threshold = pctid / 100.0

    with instream_context:
        # parse file
        for line in instream:

            # print header
            if line[0] == "@":
                print(line.rstrip())
                continue

            # get fields
            sline = line.rstrip().split('\t')
            query = sline[0]
            # code = bin(int(sline[1]))[2:].zfill(12)
            flag = int(sline[1])
            # strand = int(code[-5])
            strand = flag & 0x10
            ref = sline[2]
            cigar = sline[5]
            # cigar = HARDCLIP_RE.sub('S', cigar)
            # sline[5] = cigar
            seq = sline[9]
            qual = sline[10]

            # skip empty hits
            if ref == '*' or cigar == '*':
                continue

            # make sure read is mapped
            if flag & 0x4:
                raise ValueError("Read is unmapped.")
                # quit('ERROR 1: SAM file format')

            # calculate edit distance, total length
            hbeg, alen, hend, tlen = decode_cigar(cigar)
            mismatch = int(re.search('[NX]M:i:(\d+)', line).group(1))
            match = alen - mismatch

            # handle SAM asterisks
            if seq == '*' and qual == '*':
                # use last seq/qual
                if cquery != query:
                    raise ValueError(f"{cquery=} != {query=} in secondary alignment.")
                    # quit('ERROR 2: SAM file format')
            else:
                # update seq/qual
                cquery = query
                cseq = seq
                cqual = qual
                cstrand = strand

            # filter by percent identity / min length
            if match / tlen < id_threshold or tlen < minlen:
                continue

            # always set the seq/qual columns
            if strand == cstrand:
                sline[9] = cseq
                sline[10] = cqual
            else:
                sline[9] = reverse_complement(cseq)
                sline[10] = cqual[::-1]

            # ensure that the cigar matches the sequence
            # if tlen != (len(cseq) - hbeg - hend):
            if tlen != len(cseq):
                raise ValueError(f"Calculated {tlen=} does not match sequence length {len(cseq)} (with {cigar=}).")
                # quit('ERROR 3: SAM file format')

            # finally, print quality filtered line
            print('\t'.join(sline))


if __name__ == "__main__":
    main()
