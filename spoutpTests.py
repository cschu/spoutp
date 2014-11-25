#!/usr/bin/env python
'''
spoutp unit tests
'''

import re
import os
import sys
import subprocess as SP
from itertools import chain

import unittest

import spoutp

__author__ = 'Christian Schudoma'
__copyright__ = 'Copyright 2014, Christian Schudoma'
__license__ = 'MIT'
__version__ = '0.1a'
__maintainer__ = 'Christian Schudoma'
__email__ = 'christian.schudoma@tsl.ac.uk'

SPOUTPDIR = './'

TESTOUTPUT = open(os.path.join(SPOUTPDIR, 'test', 'multi.fa.out.short'))


class spoutpTests(unittest.TestCase):
    def test_processMultiSignalP3Output(self):
        self.failUnless(len(spoutp.processMultiSignalP3Output(TESTOUTPUT)) == 3)
        
    pass


def main(argv):
    unittest.main()
