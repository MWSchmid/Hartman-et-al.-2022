# -*- coding: utf-8 -*-
"""
usage:
python trim16SmanualBothSides.py maxMM maxOffset clipPrimerAsWell forwardSeq reverseSeq fastq(.gz)

IT IS RECOMMENDED TO REMOVE THE PRIMERS!

maxMM: maximal number of primer mismatches
maxOffset: maximal number of missing primer bases at the beginning/end
clipPrimerAsWell: everything up to the primer will be clipped. If true, then the primer as well
forwardSeq/reverseSeq: forward and reverse primer. Both will be searched on both ends:

First, forward at beginning and reverse (revCompl) at the end.
Second, reverse at beginning and forward (revCompl) at the end

note: 
* Primers must be at the end/beginning of the sequence. up to maxOffset missing bases at the beginning/end will be tolerated.
"""

import gzip
import sys

class basicRead(object):
    """A read with the headers and sequences"""
    def __init__(self, seqHead, seq, qualHead, qual):
        self.seqHead = seqHead
        self.seq = seq
        self.qualHead = qualHead
        self.qual = qual
        self.name = self.getName()
        self.matchTable = {"A":set(["A"]),
                           "C":set(["C"]),
                           "G":set(["G"]),
                           "T":set(["T"]),
                           "R":set(["A", "G"]),
                           "Y":set(["C", "T"]),
                           "S":set(["G", "C"]),
                           "W":set(["A", "T"]),
                           "K":set(["G", "T"]),
                           "M":set(["A", "C"]),
                           "B":set(["C", "G", "T"]),
                           "D":set(["A", "G", "T"]),
                           "H":set(["A", "C", "T"]),
                           "V":set(["A", "C", "G"]),
                           "N":set(["A", "C", "G", "T"])}

    def __str__(self):
        return '\n'.join([self.seqHead, self.seq, self.qualHead, self.qual])

    def getName(self):
        return self.seqHead
    
    def getSize(self):
        return len(self.seq)
    
    def chopOneAtStart(self):
        self.seq = self.seq[1:]
        self.qual = self.qual[1:]
    
    def chopOneAtEnd(self):
        self.seq = self.seq[:-1]
        self.qual = self.qual[:-1]
    
    def matchesAtStart(self, primerSeq, maxMM, doClip, toCheck):
        querySeq = list(primerSeq)
        searchSeq = list(self.seq[0:len(primerSeq)])
        if querySeq == searchSeq:
            if doClip:
                self.seq = self.seq[len(primerSeq):]
                self.qual = self.qual[len(primerSeq):]
            return True
        numMM = 0
        for (q, s) in zip(querySeq, searchSeq):
            if s not in self.matchTable[q]:
                numMM += 1
        out = numMM <= maxMM
        while (not out) and (toCheck > 0):
            toCheck -= 1
            self.chopOneAtStart()
            searchSeq = list(self.seq[0:len(primerSeq)])
            numMM = 0
            for (q, s) in zip(querySeq, searchSeq):
                if s not in self.matchTable[q]:
                    numMM += 1
            out = numMM <= maxMM
        if out and doClip:
            self.seq = self.seq[len(primerSeq):]
            self.qual = self.qual[len(primerSeq):]
        return out
        
    def matchesAtEnd(self, primerSeq, maxMM, doClip, toCheck):
        querySeq = list(primerSeq)
        searchSeq = list(self.seq[len(self.seq)-len(primerSeq):])
        if querySeq == searchSeq:
            if doClip:
                self.seq = self.seq[:len(self.seq)-len(primerSeq)]
                self.qual = self.qual[:len(self.qual)-len(primerSeq)]
            return True
        numMM = 0
        for (q, s) in zip(querySeq, searchSeq):
            if s not in self.matchTable[q]:
                numMM += 1
        out = numMM <= maxMM
        while (not out) and (toCheck > 0):
            toCheck -= 1
            self.chopOneAtEnd()
            searchSeq = list(self.seq[len(self.seq)-len(primerSeq):])
            numMM = 0
            for (q, s) in zip(querySeq, searchSeq):
                if s not in self.matchTable[q]:
                    numMM += 1
            out = numMM <= maxMM
        if out and doClip:
            self.seq = self.seq[:len(self.seq)-len(primerSeq)]
            self.qual = self.qual[:len(self.qual)-len(primerSeq)]
        return out

class illuminaRead(basicRead):
    """An illumina read with the headers and sequences"""
    def getName(self):
        return self.seqHead.split(' ')[0]


class solidRead(basicRead):
    """A solid read with the headers and sequences"""
    def getName(self):
        return '_'.join(self.seqHead.split('_')[:-1])


def myopen(fileName, mode="r"):
    """open either a regular or a compressed file"""
    if fileName.endswith(".gz"):
        return gzip.open(fileName, mode=mode)
    else:
        return open(fileName, mode=mode)


def fastqIter(infileName, readObject):
    """file iterator (returns readObjects - e.g. solidRead or illuminaRead)"""
    with myopen(infileName) as infile:
        while True:
            seqHead = infile.readline().strip()
            if not seqHead:
                break
            seq = infile.readline().strip()
            qualHead = infile.readline().strip()
            qual = infile.readline().strip()
            yield readObject(seqHead, seq, qualHead, qual)


def revcompl(seq):
    trans = {"A":"T",
             "T":"A",
             "C":"G",
             "G":"C",
             "R":"Y",
             "Y":"R",
             "K":"M",
             "M":"K",
             "S":"W",
             "W":"S",
             "B":"V",
             "D":"H",
             "H":"D",
             "V":"B",
             "N":"N"}
    revcompseq = ''.join([trans[x] for x in seq[::-1]])
    return revcompseq

if __name__ == "__main__":
    # check input
    try:
        maxMM = int(sys.argv[1])
        maxOffset = int(sys.argv[2])
        clipPrimerAsWell = sys.argv[3] == "true"
        forwardPrimer = sys.argv[4]
        reversePrimer = revcompl(sys.argv[5])
        fastqInfileName = sys.argv[6]
    except:
        print __doc__
        sys.exit(1)
    ##
    readObject = illuminaRead
    counter = 0
    hasStartMatch = 0
    hasEndMatch = 0
    hasBothMatch = 0
    totSize = 0
    ifiter = fastqIter(fastqInfileName, readObject)
    curRead = ifiter.next()
    while True:
        hasAmatch = False
        startMatch = curRead.matchesAtStart(forwardPrimer, maxMM, clipPrimerAsWell, maxOffset)
        endMatch = curRead.matchesAtEnd(reversePrimer, maxMM, clipPrimerAsWell, maxOffset)
        if startMatch:
            #hasAmatch = True
            hasStartMatch += 1
        if endMatch:
            #hasAmatch = True
            hasEndMatch += 1
        if startMatch and endMatch:
            hasAmatch = True
            hasBothMatch += 1
        if hasAmatch:
            totSize += curRead.getSize()
            print curRead
        try:
            curRead = ifiter.next()
            counter += 1
        except:
            break
        if (counter % 1000000) == 0:
            print >> sys.stderr, "processed %d million reads" % (counter)
    ## close files
    print >> sys.stderr, "%d reads in total" % counter
    print >> sys.stderr, "%d (%d / %d), %3.2f" % (hasBothMatch, hasStartMatch, hasEndMatch, float(totSize)/hasBothMatch)
