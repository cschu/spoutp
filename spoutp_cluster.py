#!/usr/bin/env python

import re
import os
import glob
import sys
import subprocess as SP
import time
import datetime

import argparse
import tempfile
import shutil

import biolib as CSBio
from spoutp import PRED_HEADER, SCORE_HEADER


# SPOUTP = 'spoutp.py'

def main(argv):
    # bsub -q TSL-Test128 "spoutp_cluster input_file num_threads<max 10>"
    # python $SPOUTP_CLUSTER $1 $2 $DEFAULT_QUEUE $PATH_TO_SPOUTP
    descr = ''
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument('--input', help='Multi-Fasta file containing CDSs.')
    parser.add_argument('--output', help='Prefix for output files.')
    parser.add_argument('--threads', help='Number of threads (max. 10).', type=int, default=10)
    parser.add_argument('--queue', default='TSL-Test128', help='Which queue to use (default: TSL-Test128).')
    parser.add_argument('--path-to-spoutp', help='Where is your installation of spoutp located?', default='/usr/users/sl/schudomc/spoutp')
    parser.add_argument('--logfile', help='Name of logfile (default: spoutp_cluster.log)', default='spoutp_cluster.log')
    args = parser.parse_args()

	
    logfile = open(args.logfile, 'wb')
    

    nSeqs = len(list(CSBio.anabl_getContigsFromFASTA(args.input)))
    nJobs = min(10, int(args.threads))

    packsize = int(float(nSeqs) / nJobs + 0.5)

    # queue = args.queue
    # path_to_spoutp = args.spoutp

    # tempfile.mkdtemp([suffix=''[, prefix='tmp'[, dir=None]]])
    # tempfile.mkstemp([suffix=''[, prefix='tmp'[, dir=None[, text=False]]]])
    tmpdir = tempfile.mkdtemp(dir='/tmp')
    tmpfiles = []

    c, i = 0, 0
    for seqid, seq in CSBio.anabl_getContigsFromFASTA(args.input):
        if c > packsize or i == 0:
            try:
                fout.close()
            except:
                pass
            # fout = open(argv[0] + '.%i' % i, 'wb')
            tmpfiles.append(tempfile.mkstemp(suffix='.%i' % i, dir=tempdir))
            fout = open(tmpfiles[-1], 'wb')
            c = 0
            i += 1
        fout.write('>%s\n%s\n' % (seqid, seq))
        c += 1
    try:
        fout.close()
    except:
        pass
    

    jobs = set()
    # for i in xrange(nJobs):
    for i, fi in enumerate(tmpfiles):
        # fi = '%s.%i' % (argv[0], i)
        cmd = 'bsub -q %s "source python-2.7.4; python %s/spoutp.py %s > %s"' % (args.queue, args.path_to_spoutp, fi, '/dev/null')
        sub = SP.Popen(cmd, shell=True, stdin=SP.PIPE, 
                       stdout=SP.PIPE, stderr=SP.PIPE)
        stdout, stderr = sub.communicate() 

        jobid = re.search('Job <([0-9]+)> is submitted', stdout)
        print stdout, stderr
        try:
            jobs.add(jobid.groups()[0])
        except:
            logfile.write('Error: job with file %s could not be submitted.\n' % fi)

    start = datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y")
    running = set(jobs)
    while True:        
        if not running: break
        if os.path.exists('KILL_SPOUTP'): sys.exit(0)
        time.sleep(60)
        logfile.write('Jobs have been running since %s. %i jobs still running\n' % (start, len(jobs)))
        
        sub = SP.Popen('bjobs', shell=True, stderr=SP.PIPE, stdout=SP.PIPE)
        stdout, stderr = sub.communicate()
        
        running = set()
        for line in stdout.split('\n')[:-1]:
            line = line.strip()
            if line.startswith('JOBID'): continue
            running.add(line.split()[0])

        running = jobs.intersection(running)

    logfile.close()
        
    summary_scores = open(args.output + '.signal_scores.tsv', 'wb')
    summary_peptides = open(args.output + '.signal_peptides.tsv', 'wb')

    summary_scores.write('\t'.join(SCORE_HEADER) + '\n')
    summary_peptides.write('\t'.join(PRED_HEADER) + '\n')
    
    # for i in xrange(nJobs):
    for i, fi in enumerate(tmpfiles):
        # fi = argv[0] + '.%i' % i
        scores = open(fi + '.signal_scores.tsv').read()
        peptides = open(fi + '.signal_peptides.tsv').read()

        summary_scores.write('\n'.join([line for line in scores.split('\n')
                                        if not line.startswith('#')]))
        summary_peptides.write('\n'.join([line for line in peptides.split('\n')
                                          if not line.startswith('#')]))

    shutil.rmtree(tmpdir)

    summary_scores.close()
    summary_peptides.close()


if __name__ == '__main__': main(sys.argv[1:])
