




import csv
import pickle
import subprocess


import numpy as np
import pandas as pd
from scipy import sparse as sps



from matplotlib import pyplot, patches
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import gridspec, colors
# from pygraphviz import *
from statannot import add_stat_annotation



def readDESeq2(file, duration):
    DF = pd.read_table(file)
    new = DF.reset_index()["index"].str.split("-", expand=True)
    DF = DF.reset_index(drop=True)
    DF["chr"] = new[0].str.strip()
    DF["start"] = new[1].str.strip().astype("int64")
    DF["end"] = new[2].str.strip().astype("int64")
    DF["duration"] = duration
    return DF




dex1 = readDESeq2("/groups/lackgrp/ll_members/berkay/StarrPipe/results/coverage/GR.1h.dex.tsv", "1h")
dex4 = readDESeq2("/groups/lackgrp/ll_members/berkay/StarrPipe/results/coverage/GR.4h.dex.tsv", "4h")
dex8 = readDESeq2("/groups/lackgrp/ll_members/berkay/StarrPipe/results/coverage/GR.8h.dex.tsv", "8h")
dex12 = readDESeq2("/groups/lackgrp/ll_members/berkay/StarrPipe/results/coverage/GR.12h.dex.tsv", "12h")



lib0 = readDESeq2("/groups/lackgrp/ll_members/berkay/StarrPipe/results/coverage/GR.0h.lib.tsv", "0h")
lib1 = readDESeq2("/groups/lackgrp/ll_members/berkay/StarrPipe/results/coverage/GR.1h.lib.tsv", "1h")
lib4 = readDESeq2("/groups/lackgrp/ll_members/berkay/StarrPipe/results/coverage/GR.4h.lib.tsv", "4h")
lib8 = readDESeq2("/groups/lackgrp/ll_members/berkay/StarrPipe/results/coverage/GR.8h.lib.tsv", "8h")
lib12 = readDESeq2("/groups/lackgrp/ll_members/berkay/StarrPipe/results/coverage/GR.12h.lib.tsv", "12h")


dex = pd.concat([dex1, dex4, dex8, dex12]).reset_index(drop=True)
dex["comparison"] = "dex"
dex["name"] = dex["chr"] +"-"+ dex["start"].map(str) +"-"+ dex["end"].map(str)


lib = pd.concat([lib0, lib1, lib4, lib8, lib12]).reset_index(drop=True)
lib["comparison"] = "lib"
lib["name"] = lib["chr"] +"-"+ lib["start"].map(str) +"-"+ lib["end"].map(str)

home = "/groups/lackgrp/ll_members/berkay/STARRbegin/peaks"

bedGR = pd.read_table(f"{home}/common.GR.peaks.bed", names=["chr", "start", "end", "name", "score"])
bedNGR = pd.read_table(f"{home}/negativeControlGR.final.bed", names=["chr", "start", "end", "name"])

def posName(bed):
    return bed["chr"] +":"+ bed["start"].map(str) +"-"+ bed["end"].map(str)

bedNGR["posName"] = posName(bedNGR)
bedGR["posName"] = posName(bedGR)

dex["posName"] = posName(dex)
lib["posName"] = posName(lib)

dex.loc[dex["posName"].isin(bedNGR["posName"]), "nodeClass"] = "nGR"
dex.loc[dex["posName"].isin(bedGR["posName"]), "nodeClass"] = "GR"

lib.loc[lib["posName"].isin(bedNGR["posName"]), "nodeClass"] = "nGR"
lib.loc[lib["posName"].isin(bedGR["posName"]), "nodeClass"] = "GR"


dex["significant"] = "no"
dex.loc[dex["padj"] < 0.05, "significant"] = "yes"

################################################################################



z = [
    (lfc, dur)
    for lfc in [0.5, 0.6, 0.7, 0.8, 0.9, 1]
    for dur in ["1h", "4h", "8h", "12h"]
]

for lfc, dur in z:









z = [
    (lfc, dur)
    for lfc in [0.5, 0.6, 0.7, 0.8, 0.9, 1]
    for dur in ["1h", "4h", "8h", "12h"]
]

N = {"lfc": [], "dur": [], "con": [], "ind": [], "non": []}

for lfc, dur in z:
    N["lfc"].append(lfc)
    N["dur"].append(dur)
    print(f"lfc:{lfc}\tdur:{dur}")
    ind = dex.loc[
        (dex["log2FoldChange"] > lfc) &
        (dex["nodeClass"] == "GR") &
        (dex["padj"] < 0.05) &
        (dex["duration"] == dur),
        ["chr", "start", "end", "posName", "nodeClass"]
    ].sort_values(["chr", "start"])
    print(f"\tIND => n:{len(ind)}")
    N["ind"].append(len(ind))
    ind["nodeClass"] = "ind"
    ind.to_csv(f"/groups/lackgrp/ll_members/berkay/StarrPipe/regions/ind.{lfc}.GR.{dur}.bed",sep="\t", index=False, header=False)
    libSub = lib[~(lib["posName"].isin(ind["posName"]))]
    con = libSub.loc[
        (libSub["log2FoldChange"] >= lfc) &
        (libSub["nodeClass"] == "GR") &
        (libSub["padj"] < 0.05) &
        (libSub["duration"] == dur),
        ["chr", "start", "end", "posName", "nodeClass"]
    ]
    print(f"\tCON => n:{len(con)}")
    N["con"].append(len(con))
    con["nodeClass"] = "con"
    non = libSub.loc[
        ~(libSub["posName"].isin(list(ind["posName"]) + list(con["posName"]))) &
        (libSub["nodeClass"] == "GR") &
        (libSub["duration"] == dur),
        ["chr", "start", "end","posName", "nodeClass"]
    ]
    print(f"\tNON => n:{len(non)}")
    N["non"].append(len(non))
    non["nodeClass"] = "non"
    con.to_csv(f"/groups/lackgrp/ll_members/berkay/StarrPipe/regions/con.{lfc}.GR.{dur}.bed",sep="\t", index=False, header=False)
    non.to_csv(f"/groups/lackgrp/ll_members/berkay/StarrPipe/regions/non.{lfc}.GR.{dur}.bed",sep="\t", index=False, header=False)


nDF = pd.DataFrame(N)


nDF.pivot()
rDF = nDF.pivot(columns="dur", index="lfc")
DF = rDF.unstack().reset_index()



fig = plt.figure(figsize=[5,5])
gs = gridspec.GridSpec(ncols=1, nrows=1)
plt.subplots_adjust(wspace=0.4)

sns.lineplot(data=nDF)

fig.savefig("starrN.pdf")
