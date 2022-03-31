#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
usage:
python splitOneFastaIntoSeveral.py fastaFile outDir maxSeqPerFile
"""

import sys

infileName = sys.argv[1]
outDir = sys.argv[2]
maxSeqPerFile = int(sys.argv[3])

seqCounter = 0
outfile = open(outDir+'/split_'+str(seqCounter)+".fasta", "wb")
with open(infileName, 'rb') as infile:
    for line in infile:
        if line[0] == '>':
            seqCounter += 1
            if (seqCounter % maxSeqPerFile) == 0:
                outfile.close()
                outfile = open(outDir+'/split_'+str(seqCounter)+".fasta", "wb")
        outfile.write(line)
outfile.close()

