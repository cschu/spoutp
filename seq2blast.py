#!/usr/bin/env python

import sys

def anabl_getContigsFromFASTA(fn, truncate=None):
    """
    Returns generator object to access sequences from a multi-FASTA file.
    Originates from 'anabl' - BLAST analysing tool, hence the prefix.
    """
    head, seq = None, ''
    for line in open(fn):
        if line[0] == '>':
            if head is not None:
                yield (head, seq) if truncate is None else (head, seq[:truncate])
            head, seq = line.strip().strip('>').split()[0], ''
        else:
            seq += line.strip()
    yield (head, seq) if truncate is None else (head, seq[:truncate])



seqs = dict(list(anabl_getContigsFromFASTA(sys.argv[1])))

for row in [line.strip() for line in open(sys.argv[2])]:
    row = row.split('\t')
    # sys.stdout.write('>%s\n%s\n' % (id_, seqs.get(id_, '')))      
    sys.stdout.write('\t'.join(row + [seqs.get(row[0], '')]) + '\n')
