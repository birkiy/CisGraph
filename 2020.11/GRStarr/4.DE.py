


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

import matplotlib
import itertools
import scipy
import math
import random
import heapq
from collections import deque


from _collections import deque
import collections
import bisect as bi




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


#
# res1hD = pd.read_table("/groups/lackgrp/ll_members/berkay/StarrPipe/results/coverage/GR.1h.dex.tsv")
# new = res1hD.reset_index()["index"].str.split("-", expand=True)
# res1hD = res1hD.reset_index(drop=True)
# res1hD["chr"] = new[0].str.strip()
# res1hD["start"] = new[1].str.strip().astype("int64")
# res1hD["end"] = new[2].str.strip().astype("int64")
# res1hD["duration"] = "1h"
#
# res4hD = pd.read_table("/groups/lackgrp/ll_members/berkay/StarrPipe/results/coverage/GR.4h.dex.tsv")
# new = res4hD.reset_index()["index"].str.split("-", expand=True)
# res4hD = res4hD.reset_index(drop=True)
# res4hD["chr"] = new[0].str.strip()
# res4hD["start"] = new[1].str.strip().astype("int64")
# res4hD["end"] = new[2].str.strip().astype("int64")
# res4hD["duration"] = "4h"
#
# res8hD = pd.read_table("/groups/lackgrp/ll_members/berkay/StarrPipe/results/coverage/GR.8h.dex.tsv")
# new = res8hD.reset_index()["index"].str.split("-", expand=True)
# res8hD = res8hD.reset_index(drop=True)
# res8hD["chr"] = new[0].str.strip()
# res8hD["start"] = new[1].str.strip().astype("int64")
# res8hD["end"] = new[2].str.strip().astype("int64")
# res8hD["duration"] = "8h"
#
# res12hD = pd.read_table("/groups/lackgrp/ll_members/berkay/StarrPipe/results/coverage/GR.12h.dex.tsv")
# new = res12hD.reset_index()["index"].str.split("-", expand=True)
# res12hD = res12hD.reset_index(drop=True)
# res12hD["chr"] = new[0].str.strip()
# res12hD["start"] = new[1].str.strip().astype("int64")
# res12hD["end"] = new[2].str.strip().astype("int64")
# res12hD["duration"] = "12h"
#



# starr = pd.concat([res0h, res1h, res4h, res8h, res12h]).reset_index(drop=True)
# starr["comparison"] = "starr"
dex = pd.concat([dex1, dex4, dex8, dex12]).reset_index(drop=True)
dex["comparison"] = "dex"

dex["name"] = dex["chr"] +"-"+ dex["start"].map(str) +"-"+ dex["end"].map(str)

dex["-log(q)"] = -1 * np.log(dex["padj"])


home = "/groups/lackgrp/ll_members/berkay/STARRbegin/peaks"

bedGR = pd.read_table(f"{home}/common.GR.peaks.bed", names=["chr", "start", "end", "name", "score"])
bedNGR = pd.read_table(f"{home}/negativeControlGR.final.bed", names=["chr", "start", "end", "name"])




def posName(bed):
    return bed["chr"] +":"+ bed["start"].map(str) +"-"+ bed["end"].map(str)


bedNGR["posName"] = posName(bedNGR)

bedGR["posName"] = posName(bedGR)


dex["posName"] = posName(dex)
dex.loc[dex["posName"].isin(bedNGR["posName"]), "nodeClass"] = "nGR"
dex.loc[dex["posName"].isin(bedGR["posName"]), "nodeClass"] = "GR"



dex["significant"] = "no"
dex.loc[dex["padj"] < 0.05, "significant"] = "yes"



z = [
    (lfc, dur)
    for lfc in [0.5, 0.6, 0.7, 0.8, 0.9, 1]
    for dur in ["1h", "4h", "8h", "12h"]
]

for lfc, dur in z:
    ind = dex.loc[
        (dex["log2FoldChange"] > lfc) &
        (dex["nodeClass"] == "GR") &
        (dex["padj"] < 0.05) &
        (dex["duration"] == dur),
        ["chr", "start", "end", "posName", "nodeClass"]
    ].sort_values(["chr", "start"])
    print(f"lfc:{lfc}\tdur:{dur}\tn:{len(ind)}")
    #nInd = dex.loc[(dex["log2FoldChange"] > 1) & (dex["nodeClass"] == "nGR") & (dex["padj"] < 0.05) & (dex["duration"] == "8h"), ["chr", "start", "end", "nodeClass"]]
    ind["nodeClass"] = "ind"
    ind.to_csv(f"/groups/lackgrp/ll_members/berkay/StarrPipe/regions/ind.{lfc}.GR.{dur}.bed",sep="\t", index=False, header=False)



























dex.loc[
    (dex["significant"] == "yes") &
    (dex["log2FoldChange"] > 0.5) &
    (dex["nodeClass"] == "GR"), "significant"] = "0.5GR"


dex.loc[
    (dex["significant"] == "yes") &
    (dex["log2FoldChange"] > 0.5) &
    (dex["nodeClass"] == "nGR"), "significant"] = "0.5nGR"

lfc05nGR = dex.loc[dex["significant"].isin(["0.5nGR"]), ["chr", "start", "end"]].drop_duplicates()
lfc05GR = dex.loc[dex["significant"].isin(["0.5GR"]), ["chr", "start", "end"]].drop_duplicates()


dex.loc[
    (dex["significant"] == "0.5GR") &
    (dex["log2FoldChange"] > 0.6) &
    (dex["nodeClass"] == "GR"), "significant"] = "0.6GR"


dex.loc[
    (dex["significant"] == "0.5nGR") &
    (dex["log2FoldChange"] > 0.6) &
    (dex["nodeClass"] == "nGR"), "significant"] = "0.6nGR"

lfc06nGR = dex.loc[dex["significant"].isin(["0.6nGR"]), ["chr", "start", "end"]].drop_duplicates()
lfc06GR = dex.loc[dex["significant"].isin(["0.6GR"]), ["chr", "start", "end"]].drop_duplicates()


dex.loc[
    (dex["significant"] == "0.6GR") &
    (dex["log2FoldChange"] > 0.7) &
    (dex["nodeClass"] == "GR"), "significant"] = "0.7GR"

dex.loc[
    (dex["significant"] == "0.6nGR") &
    (dex["log2FoldChange"] > 0.7) &
    (dex["nodeClass"] == "nGR"), "significant"] = "0.7nGR"

lfc07nGR = dex.loc[dex["significant"].isin(["0.7nGR"]), ["chr", "start", "end"]].drop_duplicates()
lfc07GR = dex.loc[dex["significant"].isin(["0.7GR"]), ["chr", "start", "end"]].drop_duplicates()


dex.loc[
    (dex["significant"] == "0.7GR") &
    (dex["log2FoldChange"] > 0.8) &
    (dex["nodeClass"] == "GR"), "significant"] = "0.8GR"

dex.loc[
    (dex["significant"] == "0.7nGR") &
    (dex["log2FoldChange"] > 0.8) &
    (dex["nodeClass"] == "nGR"), "significant"] = "0.8nGR"

lfc08nGR = dex.loc[dex["significant"].isin(["0.8nGR"]), ["chr", "start", "end"]].drop_duplicates()
lfc08GR = dex.loc[dex["significant"].isin(["0.8GR"]), ["chr", "start", "end"]].drop_duplicates()


dex.loc[
    (dex["significant"] == "0.8GR") &
    (dex["log2FoldChange"] > 0.9) &
    (dex["nodeClass"] == "GR"), "significant"] = "0.9GR"

dex.loc[
    (dex["significant"] == "0.8nGR") &
    (dex["log2FoldChange"] > 0.9) &
    (dex["nodeClass"] == "nGR"), "significant"] = "0.9nGR"

lfc09nGR = dex.loc[dex["significant"].isin(["0.9nGR"]), ["chr", "start", "end"]].drop_duplicates()
lfc09GR = dex.loc[dex["significant"].isin(["0.9GR"]), ["chr", "start", "end"]].drop_duplicates()

dex.loc[
    (dex["significant"] == "0.9GR") &
    (dex["log2FoldChange"] > 1) &
    (dex["nodeClass"] == "GR"), "significant"] = "1GR"

dex.loc[
    (dex["significant"] == "0.9nGR") &
    (dex["log2FoldChange"] > 1) &
    (dex["nodeClass"] == "nGR"), "significant"] = "1nGR"


lfc1nGR = dex.loc[dex["significant"].isin(["1nGR"]), ["chr", "start", "end"]].drop_duplicates()
lfc1GR = dex.loc[dex["significant"].isin(["1GR"]), ["chr", "start", "end"]].drop_duplicates()


import os
os.system(f"mkdir {home}/sig")
lfc05nGR.to_csv(f"{home}/sig/lfc05nGR.bed", sep="\t", index=False, header=False)
lfc05GR.to_csv(f"{home}/sig/lfc05GR.bed", sep="\t", index=False, header=False)
lfc06nGR.to_csv(f"{home}/sig/lfc06nGR.bed", sep="\t", index=False, header=False)
lfc06GR.to_csv(f"{home}/sig/lfc06GR.bed", sep="\t", index=False, header=False)
lfc07nGR.to_csv(f"{home}/sig/lfc07nGR.bed", sep="\t", index=False, header=False)
lfc07GR.to_csv(f"{home}/sig/lfc07GR.bed", sep="\t", index=False, header=False)
lfc08nGR.to_csv(f"{home}/sig/lfc08nGR.bed", sep="\t", index=False, header=False)
lfc08GR.to_csv(f"{home}/sig/lfc08GR.bed", sep="\t", index=False, header=False)
lfc09nGR.to_csv(f"{home}/sig/lfc09nGR.bed", sep="\t", index=False, header=False)
lfc09GR.to_csv(f"{home}/sig/lfc09GR.bed", sep="\t", index=False, header=False)
lfc1nGR.to_csv(f"{home}/sig/lfc1nGR.bed", sep="\t", index=False, header=False)
lfc1GR.to_csv(f"{home}/sig/lfc1GR.bed", sep="\t", index=False, header=False)


fig = plt.figure(figsize=[8,8])
gs = gridspec.GridSpec(ncols=2, nrows=2)
plt.subplots_adjust(wspace=0.3, hspace=0.3)


params =  dict(
        y="-log(q)",
        hue="significant" ,
        x="log2FoldChange",
        palette={"yes": "#660708",
                 "no":"#161A1D",
                 "1GR": "#ed6b6c",
                 "1nGR": "#856fde"}
)

for i, duration in enumerate(["1h", "4h", "8h", "12h"]):
    ax1 = fig.add_subplot(gs[i])
    sns.scatterplot(**params, data=dex[dex["duration"] == duration])
    plt.title(f"{duration} of DEX Treatment")



fig.savefig("/groups/lackgrp/ll_members/berkay/STARRbegin/results/coverage/DE.pdf")


#####################################################

ind = dex.loc[(dex["log2FoldChange"] > 1) & (dex["nodeClass"] == "GR") & (dex["padj"] < 0.05) & (dex["duration"] == "8h"), ["chr", "start", "end", "posName", "nodeClass"]]
#nInd = dex.loc[(dex["log2FoldChange"] > 1) & (dex["nodeClass"] == "nGR") & (dex["padj"] < 0.05) & (dex["duration"] == "8h"), ["chr", "start", "end", "nodeClass"]]
ind["nodeClass"] = "ind"

ind.to_csv("/groups/lackgrp/ll_members/berkay/STARRbegin/peaks/ind.GR.bed",sep="\t", index=False, header=False)


#########################33


np.mean(dex.loc[dex["nodeClass"] == "GR", "end"] - dex.loc[dex["nodeClass"] == "GR", "start"])

badBoi = dex.loc[dex["significant"].isin(["1nGR", "0.5nGR"]), ["chr", "start", "end"]]


badBoi.to_csv("/groups/lackgrp/ll_members/berkay/STARRbegin/results/peaks/significantNegativeControlGR.bed", sep="\t", index=False, header=False)










badBoi.to_csv("/groups/lackgrp/ll_members/berkay/STARRbegin/results/peaks/significantNegativeControlGR.bed", sep="\t", index=False, header=False)
