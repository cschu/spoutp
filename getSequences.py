#!/usr/bin/env python

import sys

from biolib import anabl_getContigsFromFASTA as getContigs



contigs = getContigs(sys.argv[1])

seqs = set([id_.strip().strip('>') for id_ in sys.stdin if id_])
#open('/tmp/wanted.set', 'wb').write('\n'.join(['*%s*' % seq for seq in seqs]) + '\n')
#sys.stderr.write(str(seqs) + '\n')

for id_, seq in contigs:
    # sys.stderr.write('%s %s\n' % (id_, seq))
    id_ = id_.strip().split()[0]
    if id_ in seqs:
        seqs.difference_update(set([id_]))
        sys.stdout.write('>%s\n%s\n' % (id_, seq))
    
