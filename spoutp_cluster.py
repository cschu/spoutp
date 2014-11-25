#!/usr/bin/env python

import re
import os
import glob
import sys
import subprocess as SP
import time
import datetime

import biolib as CSBio

SPOUTP = 'spoutp.py'

def main(argv):	
    logfile = open('spoutp_cluster.log', 'wb')

    nSeqs = len(list(CSBio.anabl_getContigsFromFASTA(argv[0])))
    nJobs = min(10, int(argv[1]))
    packsize = int(float(nSeqs) / nJobs + 0.5)

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
        cmd = 'bsub "source python-2.7.4; python %s %s > %s.out"' % (SPOUTP, fi, fi)
        sub = SP.Popen(cmd, shell=True, stdin=SP.PIPE, 
                       stdout=SP.PIPE, stderr=SP.PIPE)
        stdout, stderr = sub.communicate() 

        jobid = re.search('Job <([0-9]+)> is submitted', stdout)

        try:
            jobs.add(jobid.groups(0))
        except:
            logfile.write('Error: job with file %s could not be submitted.\n' % fi)

    start = datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y")
    while True:        
        if not jobs: break
        if os.path.exists('KILL_SPOUTP'): break
        time.sleep(60)
        logfile.write('Jobs have been running since %s. %i jobs still running\n' % (start, len(jobs)))
        
        sub = SP.Popen('bjobs', shell=True, stderr=SP.PIPE, stdout=SP.PIPE)
        stdout, stderr = sub.communicate()
        
        running = set()
        for line in stdout.split('\n'):
            line = line.strip()
            if line.startswith('JOBID'): continue
            running.add(line.split()[0])

        jobs = jobs - (jobs - running)
    
    logfile.close()
        
        

        





if __name__ == '__main__': main(sys.argv[1:])
