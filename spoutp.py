#!/usr/bin/env python
'''
Script to process output from SignalP-3.0.
Christian Schudoma
with appreciated support from Martin Page
'''

import re
import os
import sys
import subprocess as SP

__author__ = 'Christian Schudoma'
__copyright__ = 'Copyright 2014, Christian Schudoma'
__license__ = 'MIT'
__version__ = '0.1a'
__maintainer__ = 'Christian Schudoma'
__email__ = 'christian.schudoma@tsl.ac.uk'


VALID_START = re.compile('[ACTG]TG')

if os.uname()[1] == 'n98257':
    # local
    PATH_SIGNALP = '/home/schudomc/Downloads/SignalP-3.0/signalp-3.0/signalp'
    SRC_WRAPPER = ''
else:
    # hpc
    PATH_SIGNALP = 'signalp'
    SRC_WRAPPER = 'source signalp_wrapper-3.0_node;'
    

CODON_TABLE = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
               'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
               'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
               'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
               'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
               'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
               'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
               'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
               'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
               'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
               'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
               'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
               'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
               'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
               'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
               'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
               }

def anabl_getContigsFromFASTA(fn):
    """
    Returns generator object to access sequences from a multi-FASTA file.
    Originates from 'anabl' - BLAST analysing tool, hence the prefix.
    """
    head, seq = None, ''
    for line in open(fn):
        if line[0] == '>':
            if head is not None:
                yield (head, seq)
            head, seq = line.strip().strip('>'), ''
        else:
            seq += line.strip()
    yield (head, seq)

def getCodonsFromSequence(seq):
    """
    Returns generator object to the codons in a nucleic acid 
    sequence (or rather to trimers fron any sequence).
    Truncates the sequence to a length % 3 == 0.
    """
    it = iter(seq)
    for nt in it:
        try:
            yield nt + it.next() + it.next()
        except:
            pass
    pass

def isValidDNA(seq):
    """ Checks if a string only contains valid, non-ambiguous DNA symbols. """
    return len(re.sub('ACGT', '', seq.upper())) == 0
def isValidRNA(seq):
    """ Checks if a string only contains valid, non-ambiguous RNA symbols. """
    return len(re.sub('ACGU', '', seq.upper())) == 0

def translateCDS(seq, codonTable=CODON_TABLE):
    """ 
    Translates a CDS into a protein sequence.
    Only takes CDS with valid start codons (NUG, with BUG => AUG) as input.
    """
    peptide = ''
    checkStartCodon = True
    for codon in getCodonsFromSequence(seq.strip().upper().replace('U', 'T')):
        if checkStartCodon:
            checkStartCodon = False
            if not VALID_START.match(codon):
                sys.stderr.write('INVALID START CODON: %s\n' % codon)
                sys.exit(1)
            codon = 'ATG'
        peptide += codonTable.get(codon, 'X')
    return peptide

def callSignalP3(seqid, peptide, 
                 pathToSignalP=PATH_SIGNALP, sourceWrapper=SRC_WRAPPER):
    """
    Calls SignalP-3.0 for a single peptide sequence 
    and catches the output.
    """
    cmd = '%s%s -t euk' % (sourceWrapper, pathToSignalP)
    sub = SP.Popen(cmd, shell=True, stdin=SP.PIPE, stdout=SP.PIPE, stderr=SP.PIPE)
    stdout, stderr = sub.communicate(input='>%s\n%s\n' % (seqid, peptide))
    return stdout

def processSignalP3Output(output):
    """
    Processes SignalP-3.0 output of a single sequence. 
    """
    output = iter(output.split('\n'))
    NN_output, HMM_output = [], []
    isNN, isHMM = False, False
    while True:
        try:
            line = output.next()
        except:
            break
        line = line.strip()
        if not line: 
            continue
        if line.startswith('SignalP-NN result:'):
            isNN, isHMM = True, False
            try:
                [output.next() for i in xrange(4)]
            except:
                break
            continue
        elif line.startswith('SignalP-HMM result:'):
            isNN, isHMM = False, True
            try:
                [output.next() for i in xrange(3)]
            except:
                break            
            continue
        
        if isNN:
            if line.startswith('>'):
                try:
                    NN_result = [line] + [output.next() for i in xrange(7)]
                except:
                    break
            else:
                line = line.split()
                NN_output.append((int(line[0]), 
                                  line[1]) + tuple(map(float, line[2:])))
        elif isHMM:
            if line.startswith('>'):
                try:
                    HMM_result = [line] + [output.next() for i in xrange(4)]
                except:
                    break
            else:
                line = line.split()
                HMM_output.append((int(line[0]), 
                                  line[1]) + tuple(map(float, line[2:])))
            pass
        pass
    
    return NN_output, NN_result, HMM_output, HMM_result

def parseSignalPrediction(NN_result, HMM_result):
    """
    Parse the prediction output of SignalP-3.0,
    returning both NN and HMM results.
    """
    hmm_re = re.compile('Max cleavage site probability: (?P<p>[01]\.[0-9]{3}) between pos. (?P<start>[0-9]+) and (?P<end>[0-9]+)')
    nn_re = re.compile('# Most likely cleavage site between pos. (?P<start>[0-9]+) and (?P<end>[0-9]+):')
    
    p_hmm = hmm_re.match(HMM_result[-1])
    p_nn = nn_re.match(NN_result[-1])
    
    if p_nn:
        startNN, endNN = int(p_nn.group('start')) - 1, int(p_nn.group('end')) - 1
    else:
        startNN, endNN = None, None
    if p_hmm:
        startHMM, endHMM = int(p_hmm.group('start')) - 1, int(p_hmm.group('end')) - 1
    else:
        startHMM, endHMM = None, None
    
    return startNN, endNN, startHMM, endHMM
    

def main(argv):
    firstPositive = True
    header = ['#SeqID', 'NA_Seq', 'Pep_Seq',
              'Pos_before_CleavageSite(NN)',
              'Pos_after_CleavageSite(NN)',
              'CleavageSite(NN)',
              'NA_Seq_w/o_SigP(NN)',
              # 'Pos_before_CleavageSite(HMM)',
              # 'Pos_after_CleavageSite(HMM)',
              # 'CleavageSite(HMM)',
              # 'NA_Seq_w/o_SigP(HMM)',
              ]
    outSummary = open(argv[0] + '.signal_peptides.tsv', 'wb')

    for seqid, na_seq in anabl_getContigsFromFASTA(argv[0]):
        aa_seq = translateCDS(na_seq)        
        output = callSignalP3(seqid, aa_seq)
        NNo, NNr, HMMo, HMMr = processSignalP3Output(output)
        startNN, endNN, startHMM, endHMM = parseSignalPrediction(NNr, HMMr)

        if startNN is not None and startHMM is not None:
            if firstPositive:
                firstPositive = False
                outSummary.write('\t'.join(header) + '\n')

            # 0123456789012345678901234-56
            # MVHATSPLLLLLLLSLALVAPSLSA-RK
            #                       ***-**
            csiteNN = '%s-%s' % (aa_seq[startNN - 2:endNN], 
                                 aa_seq[endNN:endNN + 1])
            csiteHMM = '%s-%s' % (aa_seq[startHMM - 2:endHMM], 
                                  aa_seq[endHMM:endHMM + 1])
            row = [seqid, na_seq, aa_seq, 
                   endNN, endNN + 1, csiteNN, na_seq[(endNN + 1) * 3:],
                   # endHMM, endHMM + 1, csiteHMM, na_seq[(endHMM + 1) * 3:]
                   ]
            outSummary.write('\t'.join(map(str, row)) + '\n')

    outSummary.close()
    pass



if __name__ == '__main__': main(sys.argv[1:])
