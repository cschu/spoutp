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
from itertools import chain, izip

import unittest

import biolib as CSBio

__author__ = 'Christian Schudoma'
__copyright__ = 'Copyright 2014, Christian Schudoma'
__license__ = 'MIT'
__version__ = '0.1a'
__maintainer__ = 'Christian Schudoma'
__email__ = 'christian.schudoma@tsl.ac.uk'

PRED_HEADER = ['#SeqID', 'NA_Seq', 'Pep_Seq', 'seqlen_NA', 'seqlen_Pep',
               'Pos_before_CleavageSite(NN)',
               'Pos_after_CleavageSite(NN)',
               'CleavageSite(NN)',
               'NA_Seq_w/o_SigP(NN)',
               # 'Pos_before_CleavageSite(HMM)',
               # 'Pos_after_CleavageSite(HMM)',
               # 'CleavageSite(HMM)',
               # 'NA_Seq_w/o_SigP(HMM)',
               ]
SCORE_HEADER = ['#SeqID',
               'NN_Cmax_score', 'NN_Cmax_pos', 'NN_Cmax_pred',
               'NN_Ymax_score', 'NN_Ymax_pos', 'NN_Ymax_pred',
               'NN_Smax_score', 'NN_Smax_pos', 'NN_Smax_pred',
               'NN_Smean_score', 'NN_Smean_pos', 'NN_Smean_pred',
               'NN_D_score', 'NN_D_pos', 'NN_D_pred',
               'HMM_type',
               'HMM_Cmax_score', 'HMM_Cmax_pos', 'HMM_Cmax_pred',
               'HMM_Sprob_score', 'HMM_Sprob_pred'
               ]
    



if os.uname()[1] == 'n98257':
    # local
    PATH_SIGNALP = '/home/schudomc/Downloads/SignalP-3.0/signalp-3.0/signalp'
    SRC_WRAPPER = ''
else:
    # hpc
    PATH_SIGNALP = 'signalp'
    SRC_WRAPPER = 'source signalp_wrapper-3.0_node;'
    
def callMultiSignalP3(inputSequences,
                      pathToSignalP=PATH_SIGNALP, sourceWrapper=SRC_WRAPPER):
    """ 
    Calls SignalP-3.0 for a single or multiple peptide sequences
    and catches the output.
    """
    cmd = '%s%s -t euk' % (sourceWrapper, pathToSignalP)
    sub = SP.Popen(cmd, shell=True, 
                   stdin=SP.PIPE, stdout=SP.PIPE, stderr=SP.PIPE)
    sequenceData = ['>%s\n%s\n' % (seqid, peptide) 
                    for seqid, peptide in inputSequences]
    stdout, stderr = sub.communicate(input=''.join(sequenceData))
    return stdout, stderr


def processMultiSignalP3Output(output):
    """
    Processes SignalP-3.0 output of multiple sequences.
    Should also properly process single sequences, therefore 
    will replace old single sequence routine in the long run.
    Returns only the result summaries, not the scoring
    for the individual amino acids.
    Returns iterator to (seqid, NN, HMM) tuples.
    """
    output = iter(output.split('\n'))

    seqIDs, NN_gathered, HMM_gathered = [], [], []
    for line in output:
        if line.strip().startswith('-----'):
            seqIDs.append(output.next().strip().strip('>'))
        if line.strip().startswith('>'):
            currLine = line
            nextLine = output.next()
            gathered = [currLine, nextLine]
            if nextLine.strip().startswith('# Measure'):
                gathered.extend([output.next() for line in xrange(6)])
                NN_gathered.append(gathered)
            elif nextLine.strip().startswith('Prediction'):
                gathered.extend([output.next() for line in xrange(3)])
                HMM_gathered.append(gathered)
            else:
                # ignore this sequence identifier, it is outside
                # of an NN or HMM block
                pass
    # print NN_gathered, HMM_gathered
    # for x,y in izip(NN_gathered, HMM_gathered): print x, y
    return izip(seqIDs, NN_gathered, HMM_gathered)


def parseSignalPrediction(NN_result, HMM_result):
    """
    Parse the prediction output of SignalP-3.0,
    returning both NN and HMM results.
    """
    hmm_re = re.compile('Max cleavage site probability: (?P<p>[01]\.[0-9]{3}) between pos. (?P<start>[0-9]+) and (?P<end>[0-9]+)')
    nn_re = re.compile('# Most likely cleavage site between pos. (?P<start>[0-9]+) and (?P<end>[0-9]+):')
    # print 'XXXXXXXX' 
    # print NN_result
    # print 'XXXXXXXX' 
   
    # p_hmm = hmm_re.match(HMM_result)
    p_nn = nn_re.match(NN_result)

    # print 'YYYYYYYYYYY'
   
    
    if p_nn:
        # print 'ZZZZZZ', p_nn.group()
        startNN, endNN = int(p_nn.group('start')) - 1, int(p_nn.group('end')) - 1
        # print startNN, endNN, '*******'

    else:
        startNN, endNN = None, None
    #if p_hmm:
    #    startHMM, endHMM = int(p_hmm.group('start')) - 1, int(p_hmm.group('end')) - 1
    #else:
    #    startHMM, endHMM = None, None
    startHMM, endHMM = None, None
   
    
    return startNN, endNN, startHMM, endHMM


['>gi|71896386|ref|NM_0  length = 150', 
 '# Measure  Position  Value  Cutoff  signal peptide?', 
 '  max. C    26       0.838   0.32   YES', 
 '  max. Y    26       0.858   0.33   YES', 
 '  max. S    14       0.998   0.87   YES', 
 '  mean S     1-25    0.962   0.48   YES', 
 '       D     1-25    0.910   0.43   YES', 
 '# Most likely cleavage site between pos. 25 and 26: LSA-RK']

def parseSignalScores(NN_scores, HMM_scores):
    NN_re = re.compile('(?P<pos>[0-9-]+)\s+(?P<val>[01]\.[0-9]{3})\s+[01]\.[0-9]+\s+(?P<decision>Y|N)')
    NN_results = [NN_re.search(line) for line in NN_scores]
    HMM_results = ['Q' if HMM_scores[0] == 'Prediction: Non-secretory protein' 
                   else 'S']
    HMM_results.append(float(HMM_scores[-1].split(':')[1].strip().split()[0]))
    HMM_results.append(HMM_scores[-1][HMM_scores[-1].rfind(' '):])
    HMM_results.append('Y' if HMM_results[-2] > 0.5 else 'N')
    HMM_results.append(float(HMM_scores[1].split(':')[1].strip()))
    HMM_results.append('Y' if HMM_results[-1] > 0.5 else 'N')

    return list(chain.from_iterable([res.groups() if res is not None else ('','','')
                                for res in NN_results])), HMM_results


def getCleavageSiteString(seq, pos): 
    pass


def processPrediction(seqid, NNr, HMMr, filename, seqDict, outScores=None, outSummary=None):
    try:
        startNN, endNN, startHMM, endHMM = parseSignalPrediction(NNr[-1], None)#HMMr[-1])
    except:
        # print 'Failed to parse prediction'
        return outSummary, outScores
    try:
        scoresNN, scoresHMM = parseSignalScores(NNr[2:-1], HMMr[1:])
    except:
        # print 'Failed to parse scores'
        return outSummary, outScores

    if scoresNN is not None: # and scoresHMM is not None:
        if scoresHMM is None:
            scoresHMM = [''] * 6
        if outScores is None:
            outScores = open(filename + '.signal_scores.tsv', 'wb')
            outScores.write('\t'.join(SCORE_HEADER) + '\n')
        outScores.write('\t'.join(map(str, [seqid] + scoresNN + scoresHMM)) + '\n')
    if startNN is not None: # and startHMM is not None:
        if outSummary is None:
            outSummary = open(filename + '.signal_peptides.tsv', 'wb')
            outSummary.write('\t'.join(PRED_HEADER) + '\n')
        # 0123456789012345678901234-56
        # MVHATSPLLLLLLLSLALVAPSLSA-RK
        #                       ***-**
        na_seq, pep_seq = seqDict.get(seqid, ('', ''))
        csiteNN = '%s-%s' % (pep_seq[startNN - 2:endNN],
                             pep_seq[endNN:endNN + 1])
        csiteHMM = '%s-%s' % (pep_seq[startNN - 2:endHMM],
                              pep_seq[endNN:endNN + 1]) 

        row = [seqid, na_seq, pep_seq, 
               len(na_seq), len(pep_seq),
               endNN, endNN + 1, csiteNN, na_seq[endNN * 3:], pep_seq[endNN:]
               # endHMM, endHMM + 1, csiteHMM, na_seq[endHMM  * 3:]
              ]
        outSummary.write('\t'.join(map(str, row)) + '\n')

    return outSummary, outScores

def processMultiFile(filename, maxAA=70):
    firstPositive = True
    outSummary, outScores = None, None
    open(filename + '.filename', 'wb').write(filename)
    inputSequences = CSBio.translateSequences(CSBio.anabl_getContigsFromFASTA(filename))
    seqDict = {seqID: (na_seq, aa_seq) 
               for seqID, na_seq, aa_seq in inputSequences
               if aa_seq}
    open(filename + '.peptides', 'wb').write('\n'.join(['>%s\n%s' % (item[0], item[1][1]) for item in seqDict.items()]))

    input = [(seqID, seqDict[seqID][1][:maxAA].strip('*')) for seqID in seqDict] 
    open(filename + '.sigp_input', 'wb').write('\n'.join(['>%s\n%s' % item for item in input]))
  
    output, err = callMultiSignalP3(input)
    open(filename + '.sigp_output', 'wb').write(output)
    open(filename + '.sigp_errors', 'wb').write(err)
    for seqid, NNr, HMMr in processMultiSignalP3Output(output):
        outSummary, outScores = processPrediction(seqid, NNr, HMMr, 
                                                  filename, seqDict, 
                                                  outScores=outScores, 
                                                  outSummary=outSummary)
        """
        try:
            startNN, endNN, startHMM, endHMM = parseSignalPrediction(NNr[-1], 
                                                                     HMMr[-1])
        except:
            continue
        try:
            scoresNN, scoresHMM = parseSignalScores(NNr[2:-1], HMMr[1:])
        except:
            continue

        if scoresNN is not None and scoresHMM is not None:
            if outScores is None:
                outScores = open(filename + '.signal_scores.tsv', 'wb')
                outScores.write('\t'.join(SCORE_HEADER) + '\n')
            outScores.write('\t'.join(map(str, [seqid] + scoresNN + scoresHMM)) + '\n')            

        if startNN is not None and startHMM is not None:
            if outSummary is None:
                outSummary = open(filename + '.signal_peptides.tsv', 'wb')
                outSummary.write('\t'.join(PRED_HEADER) + '\n')

            # 0123456789012345678901234-56
            # MVHATSPLLLLLLLSLALVAPSLSA-RK
            #                       ***-**
            csiteNN = '%s-%s' % (seqDict[seqid][1][startNN - 2:endNN], 
                                 seqDict[seqid][1][endNN:endNN + 1])
            csiteHMM = '%s-%s' % (seqDict[seqid][1][startHMM - 2:endHMM], 
                                  seqDict[seqid][1][endHMM:endHMM + 1])
            row = [seqid, seqDict[seqid][0], seqDict[seqid][1], len(seqDict[seqid][0]), len(seqDict[seqid][1]),
                   endNN, endNN + 1, csiteNN, seqDict[seqid][0][endNN * 3:],
                   # endHMM, endHMM + 1, csiteHMM, na_seq[(endHMM + 1) * 3:]
                   ]
            outSummary.write('\t'.join(map(str, row)) + '\n')
        """

    try:
        outScores.close()
    except:
        pass
    try:
        outSummary.close()
    except:
        pass
    pass


class SpoutpTests(unittest.TestCase):
    SPOUTPDIR = './'
    try:
        TESTMULTI = open(os.path.join(SPOUTPDIR, 'test', 'multi.fa.out.short'))
        TESTSINGLE = open(os.path.join(SPOUTPDIR, 'test', 'avr2.fa.out'))
    except:
        TESTMULTI = sys.stdin
        TESTSINGLE = sys.stdin
    def test_processMultiSignalP3Output(self):
        self.failUnless(len(list(processMultiSignalP3Output(SpoutpTests.TESTMULTI.read()))) == 3)
    def test_processMultiSignalP3Output_withSingle(self):
        self.failUnless(len(list(processMultiSignalP3Output(SpoutpTests.TESTSINGLE.read()))) == 1)
    pass



def main(argv):
    if argv[0] == '--UNITTEST':
        # unittest.main()
        suite = unittest.TestLoader().loadTestsFromTestCase(SpoutpTests)
        unittest.TextTestRunner(verbosity=2).run(suite)
    elif argv[0] == '--TEST':
        # processMultiSignalP3Output(SpoutpTests.TESTOUTPUT)
        processMultiSignalP3Output(argv[1])
        outSummary, outScores = None, None
        seqDict = {}
        filename = 'testtt'
        for seqid, NNr, HMMr in processMultiSignalP3Output(open(argv[1]).read()):
            print seqid, NNr, HMMr
            outSummary, outScores = processPrediction(seqid, NNr, HMMr, 
                                                      filename, seqDict, 
                                                      outScores=outScores, 
                                                      outSummary=outSummary)
        try:
            outSummary.close()
        except:
            pass
        try:
            outScores.close()
        except:
            pass


 


    else:
        processMultiFile(argv[0])



if __name__ == '__main__': main(sys.argv[1:])
