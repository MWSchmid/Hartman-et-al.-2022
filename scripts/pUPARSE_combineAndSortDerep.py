import sys

usage = """

python %s fastaFile

""" % sys.argv[0]

if len(sys.argv) < 2:
    sys.exit(usage)

fastaFile = sys.argv[1]

def load_dna_sequences(fastafile):
    sequences = {}
    infile = open(fastafile, 'r')
    for line in infile:
        if line[0] == '>':
            name = line[0:-1].split(" ")[0]
            if name in sequences:
                print >> sys.stderr, "found identical sequence names: " + name
                name = "COPY"+name
            sequences[name] = []
        else:
            sequences[name].append(line[0:-1])
    for name in sequences:
        sequences[name] = ''.join(sequences[name])
    infile.close()
    return sequences

sequences = load_dna_sequences(fastaFile)

# check if there are identical sequences and modify the counter
seqToNames = {}
for name, seq in sequences.items():
    try:
        seqToNames[seq].append(name)
    except KeyError:
        seqToNames[seq] = [name]

combinedSum = 0
nameCounter = 0
newNameToSeq = {}
for seq, nameList in seqToNames.items():
    nameCounter += 1
    if len(nameList) > 1:
        combinedSum += 1
        newSize = 0
        for entry in nameList:
            newSize += int(entry.split('=')[-1])
    else:
        newSize = nameList[0].split('=')[-1]
    newName = ''.join([">renamed", str(nameCounter), ";size=", str(newSize)])
    newNameToSeq[newName] = seq

nameSize = [[x, int(x.split("=")[1])] for x in newNameToSeq.keys()]
nameSize.sort(lambda x,y: -cmp(x[1], y[1]))

for name, size in nameSize:
    print name
    print newNameToSeq[name]

print >> sys.stderr, "combined %d sequences" % combinedSum
