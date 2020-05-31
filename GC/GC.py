

# Script Location: CisGraph/GC/GC.py
# Utils: Vers2.0
from Functions.Utils import *

import csv
import matplotlib.pyplot as plt
import pandas as pd
import pickle
import seaborn as sns
import scipy.stats

home = "/home/ualtintas/"


def GCcontentMApping(file):
    consDict = readFasta(file)
    consID = list(consDict.keys())
    consSeq = list(consDict.values())
    i = 0
    conGC = {}
    for seqStr in consSeq:
        conGC[consID[i]] = GC(seqStr)
        i += 1
    return conGC


conGC = GCcontentMApping(home + "ARBSs/fasta/cons-enc.fasta")

indGC = GCcontentMApping(home + "ARBSs/fasta/ind.enc.fasta")

nonGC = GCcontentMApping(home + "ARBSs/fasta/nonactive-enc.fasta")

proGC = GCcontentMApping(home + "genomeAnnotations/Fasta/tss.fasta")

# proGC = pickle.load(open("Fasta/proGC.p", "rb"))

data = {
    "non": list(nonGC.values()),
    "ind": list(indGC.values()),
    "con": list(conGC.values()),
    "pro": list(proGC.values())
}




combs = [[0, 3],
[1, 3],[2, 3], [0, 2],
[1, 2], [0, 1]]
yDec = 0.04
y1 = 1.1
ylim = [0.2,1.12]


colorPalette = ["#d4f3ef", "#abf0e9", "#63b7af", "#ee8572"]

fig = plt.figure(figsize=(6,9))
boxPlot(data=data, combs=combs, colorPalette=colorPalette, y1=y1, yDec=yDec, ylim=ylim, paired=False)
plt.ylabel("GC")


fig.savefig("GC.pdf")
