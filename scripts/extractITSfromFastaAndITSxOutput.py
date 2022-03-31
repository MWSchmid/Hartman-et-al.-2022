# -*- coding: utf-8 -*-
"""
usage:
python extractITSfromFastaAndITSxOutput.py domainToExtract positionsIn fastqIn(.gz) fastqOut(.gz)
domainToExtract: ITS1
positionsIn: ITSx.positions.txt
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

def loadITSxOutput(infileName, domainToExtract):
    out = {}
    with open(infileName, 'r') as infile:
        for line in infile:
            fields = line.strip().split('\t')
            for x in fields:
                if x[:len(domainToExtract)] == domainToExtract:
                    ses = x.split(": ")[-1]
                    out[fields[0]] = [int(k) for k in ses.split('-')]
    return out

if __name__ == "__main__":
    # check input
    try:
        domainToExtract = sys.argv[1]
        positionsInfileName = sys.argv[2]
        infileName = sys.argv[3]
        outfileName = sys.argv[4]
    except:
        print __doc__
        sys.exit(1)
    toUse = loadITSxOutput(positionsInfileName, domainToExtract)
    passed = 0
    counter = 1
    with myopen(outfileName, 'w') as outfile:
        ifiter = fastaIter(infileName)
        curEntry = ifiter.next()
        while True:
            if curEntry.name in toUse:
                passed += 1
                curEntry.seq = curEntry.seq[toUse[curEntry.name][0]-1:toUse[curEntry.name][1]]
                print >> outfile, curEntry
            try:
                curEntry = ifiter.next()
                counter += 1
            except:
                break
    print >> sys.stderr, "%d / %d entries extracted" % (passed, counter)
