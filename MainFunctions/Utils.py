

from MainFunctions.Packages import *



def readBed(file):
    bed = {}
    with open(file) as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        for row in reader:
            # row3 = row[3].split(".")[0]
            row3 = row[3]
            bed[row3] = (row[0], int(row[1]),  int(row[2]))
        return bed

def sortBed(Gbed):
    Gbed = {k: v for k,v in sorted(Gbed.items(), key=lambda kv: kv[1][1])}
    Gbed = {k: v for k, v in sorted(Gbed.items(), key=lambda kv: kv[1][0])}
    return Gbed

def readFasta(file):
    consID = []
    consSeq = []
    with open(file) as fasta:
        reader = csv.reader(fasta, delimiter='\t')
        i = 0
        for row in reader:
            if (i % 2 == 0):
                consID += [row[0][1:]]
            elif (i % 2 == 1):
                consSeq += [row[0]]
            i += 1
    consDict = dict(zip(consID,consSeq))
    return consDict


def mean(list):
    return sum(list) / len(list)

def getChr(bed, chrKey):
    bed = {k: v for k, v in bed.items() if v[0] == chrKey}
    return bed

def flat(list):
    return [val for sub in list for val in sub]

def GC(seq):
    gc = {"A": 0, "T": 0, "C": 0, "G":0, "N": 0}
    for base in seq:
        gc[base.upper()] += 1
    return (gc["G"] + gc["C"]) / len(seq)

def writeBed(bed, file):
    f = open(file, "w+")
    for region in bed:
        f.write(bed[region][0]  + "\t" + str(bed[region][1]) + "\t" + str(bed[region][2]) + "\t" + region + "\n")
    f.close()

def intersect(range1, range2):
    if (range1[0] == range2[0] and
       range2[2] - range1[1] > range1[1] - range2[1] and
       range1[2] - range2[1] > range2[2] - range1[2]):
       return True
    return False
