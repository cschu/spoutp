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

def getWorkingDirectory():
    workdir = os.path.join(os.path.expanduser('~'), 'spoutp_workdir')
    if not os.path.exists(workdir):           
        os.mkdir(workdir)                     
    return tempfile.mkdtemp(dir=workdir)

def splitData(seqs, nJobs, maxSequences):
    nSeqs = len(list(seqs))
    packsize = min(maxSequences, int(float(nSeqs) / nJobs + 0.5))

    c, i = 0, 0
    for seqid, seq in seqs: 
        if c > packsize or i == 0:
            try:
                fout.close()
            except:
                pass
            tmpfiles.append(tempfile.mkstemp(suffix='.%i' % i, dir=tmpdir))
            fout = open(tmpfiles[-1][1], 'wb')
            c = 0
            i += 1
        fout.write('>%s\n%s\n' % (seqid, seq))
        c += 1
    try:
        fout.close()
    except:
        pass


def updateQueue(queue, maxJobs):
    sub = SP.Popen('bjobs', shell=True, stderr=SP.PIPE, stdout=SP.PIPE)
    stdout, stderr = sub.communicate()

    running = set()
    for line in stdout.split('\n')[:-1]:
        line = line.strip()
        if line.startswith('JOBID'): continue
        job = line.split()[0]
        if job in queue:
            running.add(job)

    queue.difference_update(queue.difference(running))
    return maxJobs - len(running) #len(queue)

def createChunk(seqs, chunksize, workdir):
    tmpfile = tempfile.mkstemp(dir=workdir)[1]
    isEmpty, hasData = False, False
    fout = open(tmpfile, 'wb')
    for i in xrange(chunksize):
        try:
            fout.write('>%s\n%s\n' % seqs.next()) 
        except:
            isEmpty = True
            break
        hasData = True
    fout.close()
    
    return tmpfile, isEmpty, hasData

def submitJob(queue, fn, hpc_queue, path_to_spoutp, logfile, is_debug):
    if is_debug:
        cmd_args = (hpc_queue, '', path_to_spoutp, fn, '/dev/null') 
    else:
        cmd_args = (hpc_queue, '-o /dev/null', path_to_spoutp, fn, '/dev/null')  
        cmd = 'bsub -q %s %s -We 480 "source python-2.7.4; python %s/spoutp.py %s > %s"' % cmd_args
    e
    logfile.write('%s\n' % cmd)
    sub = SP.Popen(cmd, shell=True, stdin=SP.PIPE, 
                   stdout=SP.PIPE, stderr=SP.PIPE)
    stdout, stderr = sub.communicate() 

    jobid = re.search('Job <([0-9]+)> is submitted', stdout)
    try:        
        queue.add(jobid.groups()[0])
    except:
        logfile.write('Error: job with file %s could not be submitted.\n' % fn)
    pass
    
def runJobs(args, maxNT=210, maxSequences=2000, maxJobs=10):
    logfile = open(args.logfile, 'wb')

    seqs = [(sid, seq[:maxNT]) 
            for sid, seq in CSBio.anabl_getContigsFromFASTA(args.input)
            if seq.startswith('ATG')]
    print 'Done reading sequences'
    nJobs = min(10, int(args.threads))

    workdir = getWorkingDirectory()
    open('workdir.name', 'wb').write(workdir)
    nSeqs = len(seqs)                                       
    chunksize = min(maxSequences, int(float(nSeqs) / nJobs + 0.5))   
    print nSeqs, chunksize
 
    queue = set()
    tmpfiles = []
    start = datetime.datetime.now()
    start_str = start.strftime("%I:%M%p on %B %d, %Y")
    it_seqs = iter(seqs)
    while True:        
        if os.path.exists('KILL_SPOUTP'): sys.exit(0)
        freeJobs = updateQueue(queue, maxJobs)
        print 'JOBS FREE: %i' % freeJobs
        for job in xrange(freeJobs):
            tmpfile, isEmpty, hasData = createChunk(it_seqs, chunksize, workdir)
            if hasData:
                tmpfiles.append(tmpfile)
                submitJob(queue, tmpfile, args.queue, args.path_to_spoutp, logfile, args.debug)
            if isEmpty:
                break
        # break
        time.sleep(60)
        pass
    duration = (start - datetime.datetime.now())#.strftime("%I:%M%p on %B %d, %Y")
    logfile.write('Start time: %s Duration: %s' % (start_str, duration))

        
    summary_scores = open(args.output + '.signal_scores.tsv', 'wb')
    summary_peptides = open(args.output + '.signal_peptides.tsv', 'wb')

    summary_scores.write('\t'.join(SCORE_HEADER) + '\n')
    summary_peptides.write('\t'.join(PRED_HEADER) + '\n')
    
    all_scores = []
    all_peptides = []

    def getKey(line):
        key = line.split()[0].strip('comp').replace('|', '_').split('_')
        return (int(key[0]), int(key[1][1:]), int(key[2][3:]), int(key[3][3:]))

    for i, fi in enumerate(tmpfiles):
        try: 
            peptides = open(fi[1] + '.signal_peptides.tsv').read()
            scores = open(fi[1] + '.signal_scores.tsv').read()
        except:
            logfile.write('Cannot open file: %s\n' % (fi[1] + '.signal_peptides.tsv'))
            continue
        all_scores.extend(sorted([(getKey(line), line)
                                  for line in scores.split('\n')
                                  if line and not line.startswith('#')]))
        all_peptides.extend(sorted([(getKey(line), line)
                                    for line in peptides.split('\n')
                                    if line and not line.startswith('#')]))


    summary_scores.write('\n'.join([item[1] for item in all_scores]))
    summary_peptides.write('\n'.join([item[1] for item in all_peptides]))

    try:
        x = 1 #shutil.rmtree(workdir)
    except:
        try:
            x = 1 #shutil.rmtree(workdir)
        except:
            logfile.write('Error deleting tempdir %s. Please remove it manually with rm -rf %s.\n' % (workdir, workdir))
        pass
    	

    logfile.close()
    summary_scores.close()
    summary_peptides.close()
    pass



def runJobsOld(args):	
    logfile = open(args.logfile, 'wb')
    
    seqs = [(sid, seq[:210]) 
            for sid, seq in CSBio.anabl_getContigsFromFASTA(args.input)
            if seq.startswith('ATG')]
    nSeqs = len(seqs)
    nJobs = min(10, int(args.threads))

    packsize = int(float(nSeqs) / nJobs + 0.5)

    # print nSeqs, nJobs, packsize
    workdir = os.path.join(os.path.expanduser('~'), 'spoutp_workdir')
   
    if not os.path.exists(workdir):
        os.mkdir(workdir)
    tmpdir = tempfile.mkdtemp(dir=workdir)
    tmpfiles = []

    c, i = 0, 0
    for seqid, seq in seqs: 
        if c > packsize or i == 0:
            try:
                fout.close()
            except:
                pass
            tmpfiles.append(tempfile.mkstemp(suffix='.%i' % i, dir=tmpdir))
            fout = open(tmpfiles[-1][1], 'wb')
            c = 0
            i += 1
        fout.write('>%s\n%s\n' % (seqid, seq))
        c += 1
    try:
        fout.close()
    except:
        pass
    

    jobs = set()
    for i, fi in enumerate(tmpfiles):
        cmd = 'bsub -q %s -We 480 "source python-2.7.4; python %s/spoutp.py %s > %s"' % (args.queue, args.path_to_spoutp, fi[1], '/dev/null')
        logfile.write('%s\n' % cmd)
        sub = SP.Popen(cmd, shell=True, stdin=SP.PIPE, 
                       stdout=SP.PIPE, stderr=SP.PIPE)
        stdout, stderr = sub.communicate() 

        jobid = re.search('Job <([0-9]+)> is submitted', stdout)
        print stdout, stderr
        try:
            jobs.add(jobid.groups()[0])
        except:
            logfile.write('Error: job with file %s could not be submitted.\n' % fi[1])

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
        
    summary_scores = open(args.output + '.signal_scores.tsv', 'wb')
    summary_peptides = open(args.output + '.signal_peptides.tsv', 'wb')

    summary_scores.write('\t'.join(SCORE_HEADER) + '\n')
    summary_peptides.write('\t'.join(PRED_HEADER) + '\n')
    
    all_scores = []
    all_peptides = []

    def getKey(line):
        key = line.split()[0].strip('comp').replace('|', '_').split('_')
        return (int(key[0]), int(key[1][1:]), int(key[2][3:]), int(key[3][3:]))

    for i, fi in enumerate(tmpfiles):
        try: 
            peptides = open(fi[1] + '.signal_peptides.tsv').read()
            scores = open(fi[1] + '.signal_scores.tsv').read()
        except:
            logfile.write('Cannot open file: %s\n' % (fi[1] + '.signal_peptides.tsv'))
            continue
        all_scores.extend(sorted([(getKey(line), line)
                                  for line in scores.split('\n')
                                  if line and not line.startswith('#')]))
        all_peptides.extend(sorted([(getKey(line), line)
                                    for line in peptides.split('\n')
                                    if line and not line.startswith('#')]))

    summary_scores.write('\n'.join([item[1] for item in all_scores]))
    summary_peptides.write('\n'.join([item[1] for item in all_peptides]))

    try:
        shutil.rmtree(tmpdir)
    except:
        try:
            shutil.rmtree(tmpdir)
        except:
            logfile.write('Error deleting tempdir %s. Please remove it manually with rm -rf %s.\n' % (tmpdir, tmpdir))
        pass
    	

    logfile.close()
    summary_scores.close()
    summary_peptides.close()
    pass


def main(argv):
    descr = ''
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument('--input', help='Multi-Fasta file containing CDSs.')
    parser.add_argument('--output', help='Prefix for output files.')
    parser.add_argument('--threads', help='Number of threads (max. 10).', type=int, default=10)
    parser.add_argument('--queue', default='TSL-Test128', help='Which queue to use (default: TSL-Test128).')
    parser.add_argument('--path-to-spoutp', help='Where is your installation of spoutp located?', default='/usr/users/sl/schudomc/spoutp')
    parser.add_argument('--logfile', help='Name of logfile (default: spoutp_cluster.log)', default='spoutp_cluster.log')
    parser.add_argument('--debug', help='Activate LSF-messages via email.', action='store_true')
    args = parser.parse_args()

    runJobs(args)
    pass







if __name__ == '__main__': main(sys.argv[1:])
