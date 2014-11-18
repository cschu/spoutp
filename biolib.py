#!/usr/bin/env python

import sys
import re
import unittest

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

VALID_START = re.compile('[ACTG]TG')


def isValidDNA(seq):
    """ Checks if a string only contains valid, non-ambiguous DNA symbols. """
    return len(re.sub('ACGT', '', seq.upper())) == 0

def isValidRNA(seq):
    """ Checks if a string only contains valid, non-ambiguous RNA symbols. """
    return len(re.sub('ACGU', '', seq.upper())) == 0

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
    sequence (or rather to consecutive trimers from any sequence).
    If length of input sequence mod 3 != 0 (i.e., the last codon is incomplete),
    sequence will be shortened by length mod 3 bases.
    """
    it = iter(seq)
    for nt in it:
        try:
            yield nt + it.next() + it.next()
        except:
            pass
    pass

def isValidStartCodon(codon, altStarts=[VALID_START]):
    return reduce(lambda x,y: x | y, 
                  map(lambda x: x.match(codon) is not None, altStarts))

def translateCDS(seq, codonTable=CODON_TABLE, altStart=[VALID_START]):
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



class CSBioTests(unittest.TestCase):
    def test_isValidXNA(self):
        self.failIf(isValidDNA('ACGJ'))
        self.failIf(isValidRNA('ACGJ'))
    def test_translateCDS(self):
        self.failUnless(translateCDS('ATG') == 'M')        
        self.failUnless(translateCDS('CTG') == 'M')
        self.failIf(translateCDS('CTGCTG') == 'MM')
        self.failUnless(translateCDS('ATGTAG') == 'M*')
        pass
    def test_getCodonsFromSequence(self):
        self.failUnless(len(list(getCodonsFromSequence('AAA'))) == 1)
        self.failUnless(len(list(getCodonsFromSequence('AAACCC'))) == 2)
        self.failUnless(len(list(getCodonsFromSequence('AAACCCGGG'))) == 3)
        self.failUnless(len(list(getCodonsFromSequence('UUUAA'))) == 1)
        pass
    def test_codonTable(self):
        self.failUnless(CODON_TABLE['TTT'] == 'F')
        self.failUnless(CODON_TABLE['TTC'] == 'F')
        self.failUnless(CODON_TABLE['TTA'] == 'L')
        self.failUnless(CODON_TABLE['TTG'] == 'L')
        self.failUnless(CODON_TABLE['TCT'] == 'S')
        self.failUnless(CODON_TABLE['TCC'] == 'S')
        self.failUnless(CODON_TABLE['TCA'] == 'S')
        self.failUnless(CODON_TABLE['TCG'] == 'S')
        self.failUnless(CODON_TABLE['TAT'] == 'Y')
        self.failUnless(CODON_TABLE['TAC'] == 'Y')
        self.failUnless(CODON_TABLE['TAA'] == '*')
        self.failUnless(CODON_TABLE['TAG'] == '*')
        self.failUnless(CODON_TABLE['TGT'] == 'C')
        self.failUnless(CODON_TABLE['TGC'] == 'C')
        self.failUnless(CODON_TABLE['TGA'] == '*')
        self.failUnless(CODON_TABLE['TGG'] == 'W')
        self.failUnless(CODON_TABLE['CTT'] == 'L')
        self.failUnless(CODON_TABLE['CTC'] == 'L')
        self.failUnless(CODON_TABLE['CTA'] == 'L')
        self.failUnless(CODON_TABLE['CTG'] == 'L')
        self.failUnless(CODON_TABLE['CCT'] == 'P')
        self.failUnless(CODON_TABLE['CCC'] == 'P')
        self.failUnless(CODON_TABLE['CCA'] == 'P')
        self.failUnless(CODON_TABLE['CCG'] == 'P')
        self.failUnless(CODON_TABLE['CAT'] == 'H')
        self.failUnless(CODON_TABLE['CAC'] == 'H')
        self.failUnless(CODON_TABLE['CAA'] == 'Q')
        self.failUnless(CODON_TABLE['CAG'] == 'Q')
        self.failUnless(CODON_TABLE['CGT'] == 'R')
        self.failUnless(CODON_TABLE['CGC'] == 'R')
        self.failUnless(CODON_TABLE['CGA'] == 'R')
        self.failUnless(CODON_TABLE['CGG'] == 'R')
        self.failUnless(CODON_TABLE['ATT'] == 'I')
        self.failUnless(CODON_TABLE['ATC'] == 'I')
        self.failUnless(CODON_TABLE['ATA'] == 'I')
        self.failUnless(CODON_TABLE['ATG'] == 'M')
        self.failUnless(CODON_TABLE['ACT'] == 'T')
        self.failUnless(CODON_TABLE['ACC'] == 'T')
        self.failUnless(CODON_TABLE['ACA'] == 'T')
        self.failUnless(CODON_TABLE['ACG'] == 'T')
        self.failUnless(CODON_TABLE['AAT'] == 'N')
        self.failUnless(CODON_TABLE['AAC'] == 'N')
        self.failUnless(CODON_TABLE['AAA'] == 'K')
        self.failUnless(CODON_TABLE['AAG'] == 'K')
        self.failUnless(CODON_TABLE['AGT'] == 'S')
        self.failUnless(CODON_TABLE['AGC'] == 'S')
        self.failUnless(CODON_TABLE['AGA'] == 'R')
        self.failUnless(CODON_TABLE['AGG'] == 'R')
        self.failUnless(CODON_TABLE['GTT'] == 'V')
        self.failUnless(CODON_TABLE['GTC'] == 'V')
        self.failUnless(CODON_TABLE['GTA'] == 'V')
        self.failUnless(CODON_TABLE['GTG'] == 'V')
        self.failUnless(CODON_TABLE['GCT'] == 'A')
        self.failUnless(CODON_TABLE['GCC'] == 'A')
        self.failUnless(CODON_TABLE['GCA'] == 'A')
        self.failUnless(CODON_TABLE['GCG'] == 'A')
        self.failUnless(CODON_TABLE['GAT'] == 'D')
        self.failUnless(CODON_TABLE['GAC'] == 'D')
        self.failUnless(CODON_TABLE['GAA'] == 'E')
        self.failUnless(CODON_TABLE['GAG'] == 'E')
        self.failUnless(CODON_TABLE['GGT'] == 'G')
        self.failUnless(CODON_TABLE['GGC'] == 'G')
        self.failUnless(CODON_TABLE['GGA'] == 'G')
        self.failUnless(CODON_TABLE['GGG'] == 'G')
        

        
        
        


def main(argv):
    unittest.main()

if __name__ == '__main__': main(sys.argv[1:])
