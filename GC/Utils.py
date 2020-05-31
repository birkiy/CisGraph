




import csv
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import pickle
import seaborn as sns
import scipy.stats
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec
import bisect as bi
import time
from statsmodels.sandbox.stats.multicomp import multipletests


def readBed(file):
    bed = {}
    with open(file) as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        for row in reader:
            # row3 = row[3].split(".")[0]
            row3 = row[3]
            bed[row3] = (row[0], int(row[1]), int(row[2]))
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





def rangesFromUpperRange(chrKey, start, end, Beds):
    values = []
    beds = {}
    for Bed in Beds:
        bed = {key: value for key, value in Beds[Bed].items() if value[0] == chrKey}
        bed = sortBed(bed)
        beds[Bed] = bed
    bedName = list(beds.keys())
    beds = list(beds.values())
    for idx, enhancer in enumerate(beds):
        enhancerEnd = [(value[2], value[1], value[0], key) for key, value in enhancer.items()]
        left = bi.bisect_left(enhancerEnd, (start, end))
        enhancerStart = [(value[1], value[2], value[0], key) for key, value in enhancer.items()]
        right = bi.bisect_right(enhancerStart, (end, start))
        values += [
            (keys, bedName[idx])
            for keys, value in list(enhancer.items())[left:right]
            if chrKey == value[0] and
               value[2] - start > start - value[1] and
               end - value[1] > value[2] - end
        ]
    return values






def NonLinCdict(steps, hexcol_array):
    cdict = {'red': (), 'green': (), 'blue': ()}
    for s, hexcol in zip(steps, hexcol_array):
        rgb =matplotlib.colors.hex2color(hexcol)
        cdict['red'] = cdict['red'] + ((s, rgb[0], rgb[0]),)
        cdict['green'] = cdict['green'] + ((s, rgb[1], rgb[1]),)
        cdict['blue'] = cdict['blue'] + ((s, rgb[2], rgb[2]),)
    return cdict



def boxPlot(data, combs, colorPalette, yDec, y1, ylim, multiTest=False, paired=False):
    PVal =[]
    for c in combs:
        if paired:
            pVal = scipy.stats.wilcoxon(data[list(data.keys())[c[0]]], data[list(data.keys())[c[1]]])[1]
        else:
            pVal = scipy.stats.ranksums(data[list(data.keys())[c[0]]], data[list(data.keys())[c[1]]])[1]
        PVal += [pVal]
    if multiTest:
        pAdj = list(multipletests(PVal, method='bonferroni')[1])
    else:
        pAdj = PVal
    cGC = flat([data[d] for d in data])
    nodeClasses = [nodeClass for nodeClass in data.keys() for _ in range(len(data[nodeClass]))]
    df2 = {
        "nodeClass": nodeClasses,
        "GC": cGC
    }
    df2 = pd.DataFrame(df2)
    sns.boxplot(x="nodeClass", y="GC", data=df2, palette=colorPalette)
    # sns.swarmplot(x="nodeClass", y="GC", data=df2, color=".25")
    plt.xlabel("")
    plt.ylabel("")
    plt.ylim(ylim)
    for c, q in zip(combs, pAdj):
        pVal = q
        x1 = c[0]
        x2 = c[1]
        # if y1 != 0.98:
        y1 = y1 - yDec
        y2 = y1 + (yDec / 3)
        if pVal < 0.0001:
            pS = "***"
        elif pVal < 0.001:
            pS = "**"
        elif pVal < 0.01:
            pS = "*"
        elif pVal < 0.05:
            pS = "."
        else:
            pS = "ns"
        plt.plot([x1, x1, x2, x2], [y1, y2, y2, y1], linewidth=1, color='k')
        plt.text((x1 + x2) * .5, y2, pS, ha='center', va='bottom', color="k", fontsize=10)
