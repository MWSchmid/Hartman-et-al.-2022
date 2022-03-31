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

nameSize = [[x, int(x.split("=")[1])] for x in sequences.keys()]
nameSize.sort(lambda x,y: -cmp(x[1], y[1]))

for name, size in nameSize:
    print name
    print sequences[name]


