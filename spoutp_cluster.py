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

CHECK_JOBQUEUE_INTERVAL = 30 # check for free jobs every 30 seconds
MAX_LEN_NT = 210 # max use first 70 aas
MAX_JOBS = 10
MAX_SEQUENCES = 1500 # wish i could increase that, but 3k, 3.5k, 4k all result in "error running HOW"

def processOutput(summary_scores, summary_peptides, seqDict, tmpfiles, workdir, logfile):

    summary_scores.write('\t'.join(SCORE_HEADER) + '\n')
    summary_peptides.write('\t'.join(PRED_HEADER) + '\n')
    
    all_scores = []
    all_peptides = []

    def getKey(line):
        key = line.split()[0].strip('comp').replace('|', '_').split('_')
        return (int(key[0]), int(key[1][1:]), int(key[2][3:]), int(key[3][3:]))

    problematic_chunks = set()
    for i, fi in enumerate(tmpfiles):
        try:
            sigp_status = open(fi + '.sigp_output').read()
        except:
            sigp_status = None
            logfile.write('Cannot open file: %s\n' % (fi + '.sig_p.output'))
        if sigp_status is None or sigp_status.startswith('error running HOW'):
            problematic_chunks.add((fi, 'HOW_error'))
 
        try: 
            peptides = open(fi + '.signal_peptides.tsv').read()
        except:
            logfile.write('Cannot open file: %s\n' % (fi + '.signal_peptides.tsv'))
            problematic_chunks.add((fi, 'no_pred'))
            continue
        try:
            scores = open(fi + '.signal_scores.tsv').read()
        except:
            logfile.write('Cannot open file: %s\n' % (fi + '.signal_scores.tsv'))
            problematic_chunks.add((fi, 'no_scores'))
            continue
         


        all_scores.extend(sorted([(getKey(line), line)
                                  for line in scores.split('\n')
                                  if line and not line.startswith('#')]))
        all_peptides.extend(sorted([(getKey(line), line)
                                    for line in peptides.split('\n')
                                    if line and not line.startswith('#')]))

    summary_scores.write('\n'.join([item[1] for item in all_scores]))
    summary_peptides.write('\n'.join([item[1] for item in all_peptides]))

    for pchunk, status in problematic_chunks:
        logfile.write('PCHUNK: %s %s\n' % (pchunk, status))
    pass


def getWorkingDirectory():
    workdir = os.path.join(os.path.expanduser('~'), 'spoutp_workdir')
    if not os.path.exists(workdir):           
        os.mkdir(workdir)                     
    return tempfile.mkdtemp(dir=workdir)

def splitData(seqs, nJobs, maxSequences, workdir):
    nSeqs = len(list(seqs))
    packsize = min(maxSequences, int(float(nSeqs) / nJobs + 0.5))

    c = 0
    for seqid, seq in seqs: 
        if c > packsize:
            try:
                fout.close()
            except:
                pass
            tmpfiles.append(tempfile.mkstemp(dir=workdir))
            fout = open(tmpfiles[-1][1], 'wb')
            c = 0
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

    # queue.difference_update(queue.difference(running)) 
    return running, maxJobs - len(running) #len(queue)

def createChunk(seqs, chunksize, workdir):#, maxNT):
    tmpfile = tempfile.mkstemp(dir=workdir)[1]
    hasData, isEmpty = False, False
    fout = open(tmpfile, 'wb')
    for i in xrange(chunksize):
        try:
            sid, seq = seqs.next()
        except:
            isEmpty = True
            break
        hasData = True
        fout.write('>%s\n%s\n' % (sid, seq))
    fout.close() 
    if not hasData:
        tmpfile = None      
 
    return tmpfile, isEmpty

def submitJob(queue, fn, hpc_queue, path_to_spoutp, logfile, is_debug):
    if is_debug:
        cmd_args = (hpc_queue, '', path_to_spoutp, fn, '/dev/null') 
    else:
        cmd_args = (hpc_queue, '-o /dev/null', path_to_spoutp, fn, '/dev/null')  
    cmd = 'bsub -q %s %s -We 480 "source python-2.7.4; python %s/spoutp.py %s > %s"' % cmd_args
    logfile.write('%s\n' % cmd)
    logfile.flush()
    sub = SP.Popen(cmd, shell=True, stdin=SP.PIPE, 
                   stdout=SP.PIPE, stderr=SP.PIPE)
    stdout, stderr = sub.communicate() 

    jobid = re.search('Job <([0-9]+)> is submitted', stdout)
    try:        
        queue.add(jobid.groups()[0])
    except:
        logfile.write('Error: job with file %s could not be submitted.\n' % fn)
        logfile.flush()
    pass
    
def runJobs(args, maxNT=MAX_LEN_NT, maxSequences=MAX_SEQUENCES, maxJobs=MAX_JOBS):
    logfile = open(args.logfile, 'wb')

    
    seqs = [(sid.split(' ')[0], seq)#[:maxNT]) 
            for sid, seq in CSBio.anabl_getContigsFromFASTA(args.input)
            if seq.startswith('ATG')]# [:35678]
    print 'Done reading sequences'
    nJobs = min(MAX_JOBS, int(args.threads))

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
    allSeqsSubmitted = False
    
    while True:        
        open(tempfile.mkstemp(dir=workdir, suffix='.queue_state')[1], 'wb').write('QUEUE:%i\n%s\n' % (len(queue), '\n'.join(queue)))
        if os.path.exists('KILL_SPOUTP'): sys.exit(0)
        queue, freeJobs = updateQueue(queue, maxJobs)
        if not allSeqsSubmitted:
            print 'JOBS FREE: %i' % freeJobs
            for job in xrange(freeJobs):
                tmpfile, allSeqsSubmitted = createChunk(it_seqs, chunksize, workdir)
                if tmpfile is not None:
                    tmpfiles.append(tmpfile)
                    submitJob(queue, tmpfile, args.queue, args.path_to_spoutp, logfile, args.get_lsf_mails)
                if allSeqsSubmitted:
                    break
        elif not queue:  
            # job manager has to run as long as there are jobs being processed
            break
        time.sleep(CHECK_JOBQUEUE_INTERVAL)
        pass
    duration = (datetime.datetime.now() - start)
    logfile.write('Start time: %s Duration: %s\n' % (start_str, duration))
    logfile.flush()

    processOutput(open(args.output + '.signal_scores.tsv', 'wb'),
                  open(args.output + '.signal_peptides.tsv', 'wb'),
                  None, #dict(seqs), # can be removed again if this works!
                  tmpfiles,
                  workdir,
                  logfile)

    try:
        x = 1 #shutil.rmtree(workdir)
    except:
        try:
            x = 1 #shutil.rmtree(workdir)
        except:
            logfile.write('Error deleting tempdir %s. Please remove it manually with rm -rf %s.\n' % (workdir, workdir))
        pass

    logfile.close()
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
    parser.add_argument('--get-lsf-mails', help='Activate LSF-messages via email.', action='store_true')
    args = parser.parse_args()

    runJobs(args)
    pass







if __name__ == '__main__': main(sys.argv[1:])
