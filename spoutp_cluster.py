#!/usr/bin/env python

import re
import os
import glob
import sys
import subprocess as SP
import time
import datetime

import biolib as CSBio
from spoutp import PRED_HEADER, SCORE_HEADER


# SPOUTP = 'spoutp.py'

def main(argv):	
    logfile = open('spoutp_cluster.log', 'wb')

    nSeqs = len(list(CSBio.anabl_getContigsFromFASTA(argv[0])))
    nJobs = min(10, int(argv[1]))

    packsize = int(float(nSeqs) / nJobs + 0.5)

    queue = argv[2] if len(argv) > 2 else 'TSL-Test128'
    path_to_spoutp = argv[3] if len(argv) > 3 else '/usr/users/sl/schudomc/spoutp'
    

    c, i = 0, 0
    for seqid, seq in CSBio.anabl_getContigsFromFASTA(argv[0]):
        if c > packsize or i == 0:
            try:
                fout.close()
            except:
                pass
            fout = open(argv[0] + '.%i' % i, 'wb')
            c = 0
            i += 1
        fout.write('>%s\n%s\n' % (seqid, seq))
        c += 1
    try:
        fout.close()
    except:
        pass
    

    jobs = set()
    for i in xrange(nJobs):
        fi = '%s.%i' % (argv[0], i)
        cmd = 'bsub -q %s "source python-2.7.4; python %s/spoutp.py %s > %s"' % (queue, path_to_spoutp, fi, '/dev/null')
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
        
    summary_scores = open(argv[0] + '.signal_scores.tsv', 'wb')
    summary_peptides = open(argv[0] + '.signal_peptides.tsv', 'wb')

    summary_scores.write('\t'.join(SCORE_HEADER) + '\n')
    summary_peptides.write('\t'.join(PRED_HEADER) + '\n')
    
    for i in xrange(nJobs):
        scores = open(argv[0] + '.%i.signal_scores.tsv' % i)
        peptides = open(argv[0] + '.%i.signal_peptides.tsv' % i)

        summary_scores.write('\n'.join(scores.read().split('\n')[1:]))
        summary_peptides.write('\n'.join(peptides.read().split('\n')[1:]))

        scores.close()
        peptides.close()
        
        os.remove(argv[0] + '.%i.signal_scores.tsv' % i)
        os.remove(argv[0] + '.%i.signal_peptides.tsv' % i)
        os.remove(argv[0] + '.%i' % i)

    summary_scores.close()
    summary_peptides.close()


if __name__ == '__main__': main(sys.argv[1:])
