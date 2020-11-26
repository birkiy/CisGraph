

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
    conCpG = {}
    conGpC = {}
    conApT = {}
    conTpA = {}
    for seqStr in consSeq:
        conGC[consID[i]] = GC(seqStr)
        conCpG[consID[i]] = NpN(seqStr, "CG")
        conGpC[consID[i]] = NpN(seqStr, "GC")
        conApT[consID[i]] = NpN(seqStr, "AT")
        conTpA[consID[i]] = NpN(seqStr, "TA")
        i += 1
    return conGC, conCpG, conGpC, conApT, conTpA


# conGC, conCpG, conGpC= GCcontentMApping(home + "ARBSs/fasta/cons-enc.fasta")

con2GC, con2CpG, con2GpC, con2ApT, con2TpA = GCcontentMApping(home + "ARBSs/fasta/con.hg19.fasta")

# indGC, indCpG, indGpC= GCcontentMApping(home + "ARBSs/fasta/ind.enc.fasta")

ind2GC, ind2CpG, ind2GpC, ind2ApT, ind2TpA = GCcontentMApping(home + "ARBSs/fasta/ind.hg19.fasta")

# nonGC, nonCpG, nonGpC= GCcontentMApping(home + "ARBSs/fasta/nonactive-enc.fasta")

non2GC, non2CpG, non2GpC, non2ApT, non2TpA = GCcontentMApping(home + "ARBSs/fasta/non.hg19.fasta")

nARGC, nARCpG, nARGpC, nARApT, nARTpA = GCcontentMApping(home + "ARBSs/fasta/negativeControl.AR.fasta")

proGC, proCpG, proGpC, proApT, proTpA = GCcontentMApping(home + "genomeAnnotations/Fasta/tss.fasta")


nonDF = pd.DataFrame({"GC": list(non2GC.values()), "nodeClass": "non"})
indDF = pd.DataFrame({"GC": list(ind2GC.values()), "nodeClass": "ind"})
conDF = pd.DataFrame({"GC": list(con2GC.values()), "nodeClass": "con"})
nARDF = pd.DataFrame({"GC": list(nARGC.values()), "nodeClass": "nAR"})
proDF = pd.DataFrame({"GC": list(proGC.values()), "nodeClass": "pro"})

GC = pd.concat([conDF, indDF, nonDF, nARDF, proDF]).reset_index(drop=True)

nonDF = pd.DataFrame({"CpG": list(non2CpG.values()), "nodeClass": "non"})
indDF = pd.DataFrame({"CpG": list(ind2CpG.values()), "nodeClass": "ind"})
conDF = pd.DataFrame({"CpG": list(con2CpG.values()), "nodeClass": "con"})
nARDF = pd.DataFrame({"CpG": list(nARCpG.values()), "nodeClass": "nAR"})
proDF = pd.DataFrame({"CpG": list(proCpG.values()), "nodeClass": "pro"})

CpG = pd.concat([conDF, indDF, nonDF, nARDF, proDF]).reset_index(drop=True)

nonDF = pd.DataFrame({"GpC": list(non2GpC.values()), "nodeClass": "non"})
indDF = pd.DataFrame({"GpC": list(ind2GpC.values()), "nodeClass": "ind"})
conDF = pd.DataFrame({"GpC": list(con2GpC.values()), "nodeClass": "con"})
nARDF = pd.DataFrame({"GpC": list(nARGpC.values()), "nodeClass": "nAR"})
proDF = pd.DataFrame({"GpC": list(proGpC.values()), "nodeClass": "pro"})

GpC = pd.concat([conDF, indDF, nonDF, nARDF, proDF]).reset_index(drop=True)


nonDF = pd.DataFrame({"ApT": list(non2ApT.values()), "nodeClass": "non"})
indDF = pd.DataFrame({"ApT": list(ind2ApT.values()), "nodeClass": "ind"})
conDF = pd.DataFrame({"ApT": list(con2ApT.values()), "nodeClass": "con"})
nARDF = pd.DataFrame({"ApT": list(nARApT.values()), "nodeClass": "nAR"})
proDF = pd.DataFrame({"ApT": list(proApT.values()), "nodeClass": "pro"})

ApT = pd.concat([conDF, indDF, nonDF, nARDF, proDF]).reset_index(drop=True)


nonDF = pd.DataFrame({"TpA": list(non2TpA.values()), "nodeClass": "non"})
indDF = pd.DataFrame({"TpA": list(ind2TpA.values()), "nodeClass": "ind"})
conDF = pd.DataFrame({"TpA": list(con2TpA.values()), "nodeClass": "con"})
nARDF = pd.DataFrame({"TpA": list(nARTpA.values()), "nodeClass": "nAR"})
proDF = pd.DataFrame({"TpA": list(proTpA.values()), "nodeClass": "pro"})


TpA = pd.concat([conDF, indDF, nonDF, nARDF, proDF]).reset_index(drop=True)


#colorPalette = ["#63b7af", "#abf0e9", "#d4f3ef", "#f5fffd", "#ee8572"]

colorPalette = {"con": "#5A5A5A",
                "ind": "#F9746D",
                "non": "#ACACAC",
                "nAR": "#F5F5F5",
                "pro": "#000000"}


fig = plt.figure(figsize=(8, 5))
gs = gridspec.GridSpec(ncols=2, nrows=1)
plt.subplots_adjust(wspace=0.4)

ax1 = fig.add_subplot(gs[0])
ax1 = sns.boxplot(data=GC, y="GC", x="nodeClass", palette=colorPalette)
plt.ylabel("GC", fontsize=14)
plt.xlabel("Node Classes", fontsize=14)
for tick in ax1.xaxis.get_major_ticks():
    tick.label.set_fontsize(13)


for tick in ax1.yaxis.get_major_ticks():
    tick.label.set_fontsize(13)



boxPairs = [("pro", "con"),
            ("con", "ind"),
            ("con", "non"),
            ("ind", "non"),
            ("con", "nAR")
            ]
add_stat_annotation(ax1, x="nodeClass", y="GC", data=GC, test="Mann-Whitney", box_pairs=boxPairs)


ax2 = fig.add_subplot(gs[1])
ax2 = sns.boxplot(data=CpG, y="CpG", x="nodeClass", palette=colorPalette)
plt.ylabel("CpG", fontsize=14)
plt.xlabel("Node Classes", fontsize=14)
for tick in ax2.xaxis.get_major_ticks():
    tick.label.set_fontsize(13)


for tick in ax2.yaxis.get_major_ticks():
    tick.label.set_fontsize(13)



add_stat_annotation(ax2, x="nodeClass", y="CpG", data=CpG, test="Mann-Whitney", box_pairs=boxPairs)

fig.savefig(f"{figureRoot}/CG.pdf")
