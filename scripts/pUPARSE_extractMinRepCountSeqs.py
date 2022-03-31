# -*- coding: utf-8 -*-
"""
usage:
python pUPARSE_extractMinRepCountSeqs.py minRepCount fastqIn(.gz) fastqOut(.gz)
"""

import gzip
import sys

class basicFasta(object):
    """a fasta entry"""
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq

    def __str__(self):
        return '\n'.join([''.join(['>',self.name]), self.seq])

def myopen(fileName, mode="r"):
    """open either a regular or a compressed file"""
    if fileName.endswith(".gz"):
        return gzip.open(fileName, mode=mode)
    else:
        return open(fileName, mode=mode)

def fastaIter(infileName):
    """file iterator (returns basicFasta objects)"""
    with myopen(infileName) as infile:
        line = infile.readline().strip()
        while True:
            name = line[1:]
            sequences = []
            line = infile.readline().strip()
            if not line:
                break
            while line[0] != '>':
                sequences.append(line)
                line = infile.readline().strip()
            seq = ''.join(sequences)
            yield basicFasta(name, seq)

if __name__ == "__main__":
    # check input
    try:
        minRepCount = int(sys.argv[1])
        infileName = sys.argv[2]
        outfileName = sys.argv[3]
    except:
        print __doc__
        sys.exit(1)
    passed = 0
    counter = 1
    with myopen(outfileName, 'w') as outfile:
        ifiter = fastaIter(infileName)
        curEntry = ifiter.next()
        while True:
            if int(curEntry.name.split('=')[-1]) >= minRepCount:
                passed += 1
                print >> outfile, curEntry
            try:
                curEntry = ifiter.next()
                counter += 1
            except:
                break
    print >> sys.stderr, "%d / %d entries extracted" % (passed, counter)
