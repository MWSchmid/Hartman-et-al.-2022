import sys
import re
import gzip

usage = """

python %s fastaFile OTUnamesToExtract

fastaFile: the genome fasta, can be gzipped
OTUnamesToExtract: one name per line

results are printed to std out

""" % sys.argv[0]

if len(sys.argv) < 2:
    sys.exit(usage)

fastaFile = sys.argv[1]
otuFile = sys.argv[2]

def myopen(fileName, mode="r"):
    """open either a regular or a compressed file"""
    if fileName.endswith(".gz"):
        return gzip.open(fileName, mode=mode)
    else:
        return open(fileName, mode=mode)

def loadFasta(fastaFile):
    genomeDict = {}
    with myopen(fastaFile) as infile:
        for line in infile:
            if line[0] == '>':
                chromosome = line.strip().split(' ')[0][1:]
                genomeDict.update({chromosome:[]})
            else:
                genomeDict[chromosome].append(line[0:-1])
        for chromosome in genomeDict:
            sequence = ''.join(genomeDict[chromosome])
            genomeDict[chromosome] = sequence
    return genomeDict

sequences = loadFasta(fastaFile)

lineCounter = 0
with myopen(otuFile) as infile:
    written = 0
    notFound = 0
    for line in infile:
        lineCounter += 1
        name = line.strip()
        try:
            curSeq = sequences[name]
            print ''.join([">", name])
            print curSeq
            written += 1
        except KeyError:
            print >> sys.stderr, name
            notFound += 1
        
print >> sys.stderr, "%d sequences were missing, %d sequences were written" % (notFound, written)
