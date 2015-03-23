#!/bin/bash

import sys

from biolib import anabl_getReadsFromFastQ as getFQ


seqs1, seqs2 = getFQ(sys.argv[1]), None #getFQ(sys.argv[2])

nSequences = 0
readsWithN = []
totalBP = 0
for id_, seq, qual in seqs1:
    nSequences += 1
    nBP = len(seq)
    totalBP += nBP
    Ncount = seq.count('N')
    if Ncount:
        readsWithN.append((Ncount, nBP))

print sys.argv[1]
print nSequences, '#reads',
print totalBP, '#bp'
print len(readsWithN), '#reads containing N (%.5f)' % (float(len(readsWithN))/nSequences)
N_content = sum([read[0] for read in readsWithN])
print N_content, '#N (%.5f)' % (float(N_content)/totalBP)
    
