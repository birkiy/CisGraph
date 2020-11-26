

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

nonDF = pd.DataFrame({"CG": list(nonGC.values()), "nodeClass": "non"})
indDF = pd.DataFrame({"CG": list(indGC.values()), "nodeClass": "ind"})
conDF = pd.DataFrame({"CG": list(conGC.values()), "nodeClass": "con"})
proDF = pd.DataFrame({"CG": list(proGC.values()), "nodeClass": "pro"})

CG = pd.concat([conDF, indDF, nonDF, proDF]).reset_index(drop=True)


colorPalette = ["#63b7af", "#abf0e9","#d4f3ef", "#ee8572"]

fig = plt.figure(figsize=(5,6))
g = sns.boxplot(data=CG, y="CG", x="nodeClass", palette=colorPalette)
plt.ylabel("CG")
boxPairs = [("pro", "con"),
            ("pro", "ind"),
            ("pro", "non"),
            ("con", "ind"),
            ("con", "non"),
            ("ind", "non")
            ]
add_stat_annotation(g, x="nodeClass", y="CG", data=CG, test="Mann-Whitney", box_pairs=boxPairs)


fig.savefig(f"{figureRoot}/CG.pdf")




combs = [[0, 3],
[1, 3],[2, 3], [0, 2],
[1, 2], [0, 1]]
yDec = 0.04
y1 = 1.1
ylim = [0.2,1.12]




fig = plt.figure(figsize=(4,6))
boxPlot(data=data, combs=combs, colorPalette=colorPalette, y1=y1, yDec=yDec, ylim=ylim, paired=False)
plt.ylabel("GC")


fig.savefig("GC.pdf")
