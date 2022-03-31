#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
usage:
python pUPARSE_splitOneFastaIntoSeveral.py fastaFile outDir maxSeqPerFile
"""

import sys

infileName = sys.argv[1]
outDir = sys.argv[2]
maxSeqPerFile = int(sys.argv[3])

sampleToCount = {}
seqCounter = 0
outfile = open(outDir+'/temp_'+str(seqCounter)+".fasta", "wb")
with open(infileName, 'rb') as infile:
    for line in infile:
        if line[0] == '>':
            #
            name = line[1:-1]
            (sample, readNumber) = name.split('.')
            try:
                sampleToCount[sample] += 1
            except KeyError:
                sampleToCount[sample] = 1
            line = ''.join(['>', sample, '.', str(sampleToCount[sample]), '\n']) 
            #
            seqCounter += 1
            if (seqCounter % maxSeqPerFile) == 0:
                outfile.close()
                #
                for sample in sampleToCount:
                    sampleToCount[sample] = 0
                #
                outfile = open(outDir+'/temp_'+str(seqCounter)+".fasta", "wb")
        outfile.write(line)
outfile.close()

