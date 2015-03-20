#!/usr/bin/python2.7

import sys
import os

if len(sys.argv) < 2:
    print "Usage: %s <fastq.filename> [crappy]" % sys.argv[0]
    print "Example: %s bac.fastq" % sys.argv[0]
    sys.exit()

if len(sys.argv) == 2:
   os.system("bin/prepare_cr_index %s tmp.index" % sys.argv[1])
else:
   os.system("bin/prepare_cr_index %s tmp.index crappy" % sys.argv[1])

os.system("bin/test_cr_index tmp.index")
os.system("rm tmp.index")
