


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


starr = pd.read_table("/groups/lackgrp/ll_members/berkay/STARRbegin/totalCountTable.tsv")

inputs = ["'input.GSE114063.pool10.final.bam'",
          "'input.GSE114063.pool11.final.bam'",
          "'input.GSE114063.pool12.final.bam'",
          "'input.GSE114063.pool1.final.bam'",
          "'input.GSE114063.pool2.final.bam'",
          "'input.GSE114063.pool3.final.bam'",
          "'input.GSE114063.pool4.final.bam'",
          "'input.GSE114063.pool5.final.bam'",
          "'input.GSE114063.pool6.final.bam'",
          "'input.GSE114063.pool7.final.bam'",
          "'input.GSE114063.pool8.final.bam'",
          "'input.GSE114063.pool9.final.bam'"]
reps0h = ["'GR.0h.dex.rep1.final.bam'",
          "'GR.0h.dex.rep2.final.bam'",
          "'GR.0h.dex.rep3.final.bam'",
          "'GR.0h.dex.rep4.final.bam'",
          "'GR.0h.dex.rep5.final.bam'"]


reps12h = ["'GR.12h.dex.rep1.final.bam'",
          "'GR.12h.dex.rep2.final.bam'",
          "'GR.12h.dex.rep3.final.bam'",
          "'GR.12h.dex.rep4.final.bam'",
          "'GR.12h.dex.rep5.final.bam'"]

reps1h = ["'GR.1h.dex.rep1.final.bam'",
          "'GR.1h.dex.rep2.final.bam'",
          "'GR.1h.dex.rep3.final.bam'",
          "'GR.1h.dex.rep4.final.bam'",
          "'GR.1h.dex.rep5.final.bam'"]


reps4h = ["'GR.4h.dex.rep1.final.bam'",
       "'GR.4h.dex.rep2.final.bam'",
       "'GR.4h.dex.rep3.final.bam'",
       "'GR.4h.dex.rep4.final.bam'",
       "'GR.4h.dex.rep5.final.bam'"]

reps8h = ["'GR.8h.dex.rep1.final.bam'",
       "'GR.8h.dex.rep2.final.bam'",
       "'GR.8h.dex.rep3.final.bam'",
       "'GR.8h.dex.rep4.final.bam'",
       "'GR.8h.dex.rep5.final.bam'"]



starr["inputSum"] = np.sum(starr[inputs], axis=1)
starr["0hSum"] = np.sum(starr[reps0h], axis=1)
starr["1hSum"] = np.sum(starr[reps1h], axis=1)
starr["4hSum"] = np.sum(starr[reps4h], axis=1)
starr["8hSum"] = np.sum(starr[reps8h], axis=1)
starr["12hSum"] = np.sum(starr[reps12h], axis=1)







starr["Sum"] = np.sum(starr[reps0h], axis=1)
h0 = starr[["inputSum", "Sum"]]
h0["duration"] = "0h"
starr["Sum"] = np.sum(starr[reps1h], axis=1)
h1 = starr[["inputSum", "Sum"]]
h1["duration"] = "1h"
starr["Sum"] = np.sum(starr[reps4h], axis=1)
h4 = starr[["inputSum", "Sum"]]
h4["duration"] = "4h"
starr["Sum"] = np.sum(starr[reps8h], axis=1)
h8 = starr[["inputSum", "Sum"]]
h8["duration"] = "8h"
starr["Sum"] = np.sum(starr[reps12h], axis=1)
h12 = starr[["inputSum", "Sum"]]
h12["duration"] = "12h"

sumDF = pd.concat([h0,h1,h4,h8,h12])


def abline(slope, intercept):
    """Plot a line from slope and intercept"""
    axes = plt.gca()
    x_vals = np.array(axes.get_xlim())
    y_vals = intercept + slope * x_vals
    plt.plot(x_vals, y_vals, '--')


colorPalette = ["#FFDAB9","#FBC4AB","#F8AD9D","#F4978E","#F08080"]

fig = plt.figure(figsize=[16,9])
gs = gridspec.GridSpec(ncols=5, nrows=1)
plt.subplots_adjust(wspace=0.4, hspace=0.4)

params = dict(
    data=starr,
    x="inputSum",

)

# plt.title(f"DEX treatment for {condition} hours")
for i, (y, color) in enumerate(zip(["0hSum", "1hSum", "4hSum", "8hSum", "12hSum"], colorPalette)):
    ax1 = fig.add_subplot(gs[i])
    sns.regplot(**params, y=y, color=color)

fig.savefig("/groups/lackgrp/ll_members/berkay/STARRbegin/results/coverage/GR/GRinuptvsReg.pdf")



fig = plt.figure(figsize=[8,8])

# plt.title(f"DEX treatment for {condition} hours")
g = sns.lmplot(data=sumDF, x="inputSum", y="Sum", hue="duration", palette=colorPalette)

g.savefig("/groups/lackgrp/ll_members/berkay/STARRbegin/results/coverage/GR/GRinuptvsSum.pdf")

#
# res0h = pd.read_table("/groups/lackgrp/ll_members/berkay/STARRbegin/results/coverage/GR/GR.0h.deg.tsv")
# new = res0h.reset_index()["index"].str.split("-", expand=True)
# res0h = res0h.reset_index(drop=True)
# res0h["chr"] = new[0].str.strip()
# res0h["start"] = new[1].str.strip().astype("int64")
# res0h["end"] = new[2].str.strip().astype("int64")
# res0h["duration"] = "0h"
#
# res1h = pd.read_table("/groups/lackgrp/ll_members/berkay/STARRbegin/results/coverage/GR/GR.1h.deg.tsv")
# new = res1h.reset_index()["index"].str.split("-", expand=True)
# res1h = res1h.reset_index(drop=True)
# res1h["chr"] = new[0].str.strip()
# res1h["start"] = new[1].str.strip().astype("int64")
# res1h["end"] = new[2].str.strip().astype("int64")
# res1h["duration"] = "1h"
#
# res4h = pd.read_table("/groups/lackgrp/ll_members/berkay/STARRbegin/results/coverage/GR/GR.4h.deg.tsv")
# new = res4h.reset_index()["index"].str.split("-", expand=True)
# res4h = res4h.reset_index(drop=True)
# res4h["chr"] = new[0].str.strip()
# res4h["start"] = new[1].str.strip().astype("int64")
# res4h["end"] = new[2].str.strip().astype("int64")
# res4h["duration"] = "4h"
#
# res8h = pd.read_table("/groups/lackgrp/ll_members/berkay/STARRbegin/results/coverage/GR/GR.8h.deg.tsv")
# new = res8h.reset_index()["index"].str.split("-", expand=True)
# res8h = res8h.reset_index(drop=True)
# res8h["chr"] = new[0].str.strip()
# res8h["start"] = new[1].str.strip().astype("int64")
# res8h["end"] = new[2].str.strip().astype("int64")
# res8h["duration"] = "8h"
#
# res12h = pd.read_table("/groups/lackgrp/ll_members/berkay/STARRbegin/results/coverage/GR/GR.12h.deg.tsv")
# new = res12h.reset_index()["index"].str.split("-", expand=True)
# res12h = res12h.reset_index(drop=True)
# res12h["chr"] = new[0].str.strip()
# res12h["start"] = new[1].str.strip().astype("int64")
# res12h["end"] = new[2].str.strip().astype("int64")
# res12h["duration"] = "12h"
#########################

dexDE = dex[(dex["padj"] < 0.5)]

header = ["chr", "start", "end",
    "GR.0h.dex.rep1", "GR.0h.dex.rep2", "GR.0h.dex.rep3", "GR.0h.dex.rep4", "GR.0h.dex.rep5",
    "GR.12h.dex.rep1", "GR.12h.dex.rep2", "GR.12h.dex.rep3", "GR.12h.dex.rep4", "GR.12h.dex.rep5",
    "GR.1h.dex.rep1", "GR.1h.dex.rep2", "GR.1h.dex.rep3", "GR.1h.dex.rep4", "GR.1h.dex.rep5",
    "GR.4h.dex.rep1", "GR.4h.dex.rep2", "GR.4h.dex.rep3", "GR.4h.dex.rep4", "GR.4h.dex.rep5",
    "GR.8h.dex.rep1", "GR.8h.dex.rep2", "GR.8h.dex.rep3", "GR.8h.dex.rep4", "GR.8h.dex.rep5",
    "input.GSE114063.pool10", "input.GSE114063.pool11","input.GSE114063.pool12", "input.GSE114063.pool1", "input.GSE114063.pool2", "input.GSE114063.pool3", "input.GSE114063.pool4","input.GSE114063.pool5","input.GSE114063.pool6","input.GSE114063.pool7","input.GSE114063.pool8","input.GSE114063.pool9"]

CT = pd.read_table("/groups/lackgrp/ll_members/berkay/STARRbegin/totalCountTable.tsv", names=header, skiprows=1)


CT2 = pd.DataFrame()

CT2["input"] = CT[["input.GSE114063.pool10", "input.GSE114063.pool11","input.GSE114063.pool12", "input.GSE114063.pool1", "input.GSE114063.pool2", "input.GSE114063.pool3", "input.GSE114063.pool4","input.GSE114063.pool5","input.GSE114063.pool6","input.GSE114063.pool7","input.GSE114063.pool8","input.GSE114063.pool9"]].sum(axis=1)

CT2["0h"] = (CT[["GR.0h.dex.rep1", "GR.0h.dex.rep2", "GR.0h.dex.rep3", "GR.0h.dex.rep4", "GR.0h.dex.rep5"]].sum(axis=1) +1) / CT2["input"]
CT2["1h"] = (CT[["GR.1h.dex.rep1", "GR.1h.dex.rep2", "GR.1h.dex.rep3", "GR.1h.dex.rep4", "GR.1h.dex.rep5"]].sum(axis=1) +1) / CT2["input"]
CT2["4h"] = (CT[["GR.4h.dex.rep1", "GR.4h.dex.rep2", "GR.4h.dex.rep3", "GR.4h.dex.rep4", "GR.4h.dex.rep5"]].sum(axis=1) +1) / CT2["input"]
CT2["8h"] = (CT[["GR.8h.dex.rep1", "GR.8h.dex.rep2", "GR.8h.dex.rep3", "GR.8h.dex.rep4", "GR.8h.dex.rep5"]].sum(axis=1) +1) / CT2["input"]
CT2["12h"] = (CT[["GR.12h.dex.rep1", "GR.12h.dex.rep2", "GR.12h.dex.rep3", "GR.12h.dex.rep4", "GR.12h.dex.rep5"]].sum(axis=1) +1) / CT2["input"]


sum(np.log2(CT2.iloc[:,2] / CT2.iloc[:,1]))






# GR = pd.concat([starr, dex]).reset_index(drop=True)
# GR["name"] = GR["chr"] +"-"+ GR["start"].map(str) +"-"+ GR["end"].map(str)
#
#
# common = (GR.groupby("name").size() == 9).reset_index()
# common = list(common.loc[common[0] == True, "name"])
#
# GR = GR[GR["name"].isin(common)].sort_values("name")


names = GR.groupby(["duration", "comparison"]).get_group(("0h", "starr"))["name"].reset_index(drop=True)
lfc0s = GR.groupby(["duration", "comparison"]).get_group(("0h", "starr"))["log2FoldChange"].reset_index(drop=True)
q0s = GR.groupby(["duration", "comparison"]).get_group(("0h", "starr"))["padj"].reset_index(drop=True)

lfc1s = GR.groupby(["duration", "comparison"]).get_group(("1h", "starr"))["log2FoldChange"].reset_index(drop=True)
q1s = GR.groupby(["duration", "comparison"]).get_group(("1h", "starr"))["padj"].reset_index(drop=True)
lfc1d = GR.groupby(["duration", "comparison"]).get_group(("1h", "dex"))["log2FoldChange"].reset_index(drop=True)
q1d = GR.groupby(["duration", "comparison"]).get_group(("1h", "dex"))["padj"].reset_index(drop=True)

lfc4s = GR.groupby(["duration", "comparison"]).get_group(("4h", "starr"))["log2FoldChange"].reset_index(drop=True)
q4s = GR.groupby(["duration", "comparison"]).get_group(("4h", "starr"))["padj"].reset_index(drop=True)
lfc4d = GR.groupby(["duration", "comparison"]).get_group(("4h", "dex"))["log2FoldChange"].reset_index(drop=True)
q4d = GR.groupby(["duration", "comparison"]).get_group(("4h", "dex"))["padj"].reset_index(drop=True)

lfc8s = GR.groupby(["duration", "comparison"]).get_group(("8h", "starr"))["log2FoldChange"].reset_index(drop=True)
q8s = GR.groupby(["duration", "comparison"]).get_group(("8h", "starr"))["padj"].reset_index(drop=True)
lfc8d = GR.groupby(["duration", "comparison"]).get_group(("8h", "dex"))["log2FoldChange"].reset_index(drop=True)
q8d = GR.groupby(["duration", "comparison"]).get_group(("8h", "dex"))["padj"].reset_index(drop=True)

lfc12s = GR.groupby(["duration", "comparison"]).get_group(("12h", "starr"))["log2FoldChange"].reset_index(drop=True)
q12s = GR.groupby(["duration", "comparison"]).get_group(("12h", "starr"))["padj"].reset_index(drop=True)
lfc12d = GR.groupby(["duration", "comparison"]).get_group(("12h", "dex"))["log2FoldChange"].reset_index(drop=True)
q12d = GR.groupby(["duration", "comparison"]).get_group(("12h", "dex"))["padj"].reset_index(drop=True)

tGR = pd.DataFrame({"name": names, "lfc0s": lfc0s, "q0s": q0s,
    "lfc1s": lfc1s, "q1s": q1s, "lfc1d": lfc1d, "q1d": q1d,
    "lfc4s": lfc4s, "q4s": q4s, "lfc4d": lfc4d, "q4d": q4d,
    "lfc8s": lfc8s, "q8s": q8s, "lfc8d": lfc8d, "q8d": q8d,
    "lfc12s": lfc12s, "q12s": q12s, "lfc12d": lfc12d, "q12d": q12d})


tGR[((tGR["q1d"] < 0.05) & (tGR["lfc1d"] > 1))]
tGR[((tGR["q4d"] < 0.05) & (tGR["lfc4d"] > 1))]
tGR[((tGR["q8d"] < 0.05) & (tGR["lfc8d"] > 1))]
tGR[((tGR["q12d"] < 0.05) & (tGR["lfc12d"] > 1))]


ind = set(
    list(tGR.loc[(
        (tGR["q1d"] < 0.05) & (tGR["lfc1d"] > 1) &
        (tGR["lfc1s"] > 1) & (tGR["q1s"] < 0.05)
        ), "name"]) +
    list(tGR.loc[(
        (tGR["q4d"] < 0.05) & (tGR["lfc4d"] > 1) &
        (tGR["lfc4s"] > 1) & (tGR["q4s"] < 0.05)
        ), "name"]) +
    list(tGR.loc[(
        (tGR["q8d"] < 0.05) & (tGR["lfc8d"] > 1) &
        (tGR["lfc8s"] > 1) & (tGR["q8s"] < 0.05)
        ), "name"]) +
    list(tGR.loc[(
        (tGR["q12d"] < 0.05) & (tGR["lfc12d"] > 1) &
        (tGR["lfc12s"] > 1) & (tGR["q12s"] < 0.05)
        ), "name"])
    )

con = set(
    list(tGR.loc[(
        (tGR["q0s"] < 0.05) & (tGR["lfc0s"] > 1) &
        (tGR["q1s"] < 0.05) & (tGR["lfc1s"] > 1) &
        (tGR["lfc1d"] < 1) & (tGR["lfc1d"] > -1)
        ), "name"])
    +
    list(tGR.loc[(
        (tGR["q0s"] < 0.05) & (tGR["lfc0s"] > 1) &
        (tGR["q4s"] < 0.05) & (tGR["lfc4s"] > 1) &
        (tGR["lfc4d"] < 1) & (tGR["lfc4d"] > -1)
        ), "name"])
    +
    list(tGR.loc[(
        (tGR["q0s"] < 0.05) & (tGR["lfc0s"] > 1) &
        (tGR["q8s"] < 0.05) & (tGR["lfc8s"] > 1) &
        (tGR["lfc8d"] < 1) & (tGR["lfc8d"] > -1)
        ), "name"])
    +
    list(tGR.loc[(
        (tGR["q0s"] < 0.05) & (tGR["lfc0s"] > 1) &
        (tGR["q12s"] < 0.05) & (tGR["lfc12s"] > 1) &
        (tGR["lfc12d"] < 1) & (tGR["lfc12d"] > -1)
        ), "name"])
    )






non = set(
    list(tGR.loc[(
        (tGR["q0s"] < 0.05) & (tGR["lfc0s"] < 1) & (tGR["lfc0s"] > -1) &
        (tGR["q1s"] < 0.05) & (tGR["lfc1s"] < 1) & (tGR["lfc1s"] > -1) &
        (tGR["lfc1d"] < 1) & (tGR["lfc1d"] > -1)
        ), "name"])
    +
    list(tGR.loc[(
        (tGR["q0s"] < 0.05) & (tGR["lfc0s"] < 1) & (tGR["lfc0s"] > -1) &
        (tGR["q4s"] < 0.05) & (tGR["lfc4s"] < 1) & (tGR["lfc4s"] > -1) &
        (tGR["lfc4d"] < 1) & (tGR["lfc4d"] > -1)
        ), "name"])
    +
    list(tGR.loc[(
        (tGR["q0s"] < 0.05) & (tGR["lfc0s"] < 1) & (tGR["lfc0s"] > -1) &
        (tGR["q8s"] < 0.05) & (tGR["lfc8s"] < 1) & (tGR["lfc8s"] > -1) &
        (tGR["lfc8d"] < 1) & (tGR["lfc8d"] > -1)
        ), "name"])
    +
    list(tGR.loc[(
        (tGR["q0s"] < 0.05) & (tGR["lfc0s"] < 1) & (tGR["lfc0s"] > -1) &
        (tGR["q12s"] < 0.05) & (tGR["lfc12s"] < 1) & (tGR["lfc12s"] > -1) &
        (tGR["lfc12d"] < 1) & (tGR["lfc12d"] > -1)
        ), "name"])
    )


["ind" for i in range(len(ind))] +
["con" for i in range(len(con))] +
["non" for i in range(len(non))]

sigEnh = pd.DataFrame({
    "name": list(ind) + list(con) + list(non),
    "nodeClass":
        ["ind" for i in range(len(ind))] +
        ["con" for i in range(len(con))] +
        ["non" for i in range(len(non))]
    })



bed = sigEnh["name"].str.split("-", expand=True)
bed[3] = bed[0].str.split("chr", expand=True)[1].map(str) +":"+ (bed[1].astype("int64") + (bed[2].astype("int64") - bed[1].astype("int64"))//2).map(str)

bed.to_csv("/groups/lackgrp/ll_members/berkay/STARRbegin/results/coverage/GR/GRde.bed", header=False, index=False, sep="\t")

# ind
indbed = pd.DataFrame({
    "name": list(ind) ,
    "nodeClass":
        ["ind" for i in range(len(ind))]
    })

indbed = indbed["name"].str.split("-", expand=True)
indbed[3] = indbed[0].str.split("chr", expand=True)[1].map(str) +":"+ (indbed[1].astype("int64") + (indbed[2].astype("int64") - indbed[1].astype("int64"))//2).map(str)

indbed.to_csv("/groups/lackgrp/ll_members/berkay/STARRbegin/results/coverage/GR/GRdeInd.bed", header=False, index=False, sep="\t")

# con
conbed = pd.DataFrame({
    "name": list(con) ,
    "nodeClass":
        ["ind" for i in range(len(con))]
    })

conbed = conbed["name"].str.split("-", expand=True)
conbed[3] = conbed[0].str.split("chr", expand=True)[1].map(str) +":"+ (conbed[1].astype("int64") + (conbed[2].astype("int64") - conbed[1].astype("int64"))//2).map(str)

conbed.to_csv("/groups/lackgrp/ll_members/berkay/STARRbegin/results/coverage/GR/GRdeCon.bed", header=False, index=False, sep="\t")


# non
nonbed = pd.DataFrame({
    "name": list(non) ,
    "nodeClass":
        ["ind" for i in range(len(non))]
    })

nonbed = nonbed["name"].str.split("-", expand=True)
nonbed[3] = nonbed[0].str.split("chr", expand=True)[1].map(str) +":"+ (nonbed[1].astype("int64") + (nonbed[2].astype("int64") - nonbed[1].astype("int64"))//2).map(str)

nonbed.to_csv("/groups/lackgrp/ll_members/berkay/STARRbegin/results/coverage/GR/GRdeNon.bed", header=False, index=False, sep="\t")





tGR["nodeClass"] = "N.S."
tGR.loc[tGR["name"].isin(non), "nodeClass"] = "non"
tGR.loc[tGR["name"].isin(con), "nodeClass"] = "con"
tGR.loc[tGR["name"].isin(ind), "nodeClass"] = "ind"






fig = plt.figure(figsize=[8,8])
gs = gridspec.GridSpec(ncols=2, nrows=2)
plt.subplots_adjust(wspace=0.4, hspace=0.4)

params = dict(
    data=tGR,
    hue="nodeClass",
    s=8
)

X = ["lfc1s", "lfc4s", "lfc8s", "lfc12s"]
Y = ["lfc1d", "lfc4d", "lfc8d", "lfc12d"]
conditions = [1, 4, 8, 12]
for i,x,y,condition in zip(range(4),X,Y,conditions):
    ax1 = fig.add_subplot(gs[i])
    sns.scatterplot(**params, x=x, y=y)
    plt.title(f"DEX treatment for {condition} hours")

fig.savefig("/groups/lackgrp/ll_members/berkay/STARRbegin/results/coverage/GR/GRscattNodeClass.pdf")



labels = ["ind", "con", "non"]
sizes = [
    [3, 45, 1344],
    [19, 48, 1234],
    [55, 49, 1152],
    [59, 42, 1110]
    ]

fig = plt.figure(figsize=[8,8])
gs = gridspec.GridSpec(ncols=2, nrows=2)
plt.subplots_adjust(wspace=0.4)

conditions = [1, 4, 8, 12]
for i, condition, size in zip(range(4),conditions, sizes):
    ax1 = fig.add_subplot(gs[i])
    plt.pie(size, labels=labels, autopct='%1.0f%%')
    plt.title(f"DEX treatment for {condition} hours\n(n={sum(size)})")

fig.savefig(f"/groups/lackgrp/ll_members/berkay/STARRbegin/results/coverage/GR/GRenhPie.pdf")









tGR[(
    ((tGR["q0s"] < 0.1) & (tGR["lfc0s"] > 1)) &
    ((tGR["q4s"] < 0.1) & (tGR["lfc4s"] > 1)) &
    ((tGR["lfc4d"] < 1) & (tGR["lfc4d"] > -1))
    )]


tGR[(
    ((tGR["q0s"] < 0.05) & (tGR["lfc0s"] < 1)) &
    ((tGR["q4s"] < 0.05) & (tGR["lfc4s"] < 1)) &
    ((tGR["lfc4d"] < 1) & (tGR["lfc4d"] > -1))
)]




t




















&
(tGR["lfc12d"] < 1) & (tGR["lfc12d"] > -1)



(tGR["q0s"] < 0.05) &
(tGR["q0s"] < 0.05) &
(tGR["q0s"] < 0.05) &
(tGR["q0s"] < 0.05) &


ind = set(GR.loc[
    (GR["padj"] < 0.05) &
    ((GR["log2FoldChange"] > 1) & (GR["comparison"] == "dex")),
    "name"
    ])

con = set(GR.loc[
    (GR["padj"] < 0.05) &
    ((GR["log2FoldChange"] > 1) & (GR["comparison"] == "starr")) &
    ((GR["log2FoldChange"] < 1) & (GR["comparison"] == "dex")),
    "name"
    ])


con = set(GR.loc[
    (GR["padj"] < 0.05) &
    ((GR["log2FoldChange"] > 1) & (GR["comparison"] == "starr")),
    "name"
    ])

GRm = GR.set_index(["duration", "name"])

GRm = GRm.sort_index()


bm0h = list(GRm.loc[GRm.index.get_level_values(0) == "0h", "baseMean"])

GRm.loc[GRm.index.get_level_values(0) == "1h", "baseMean"] / bm0h

GRm.loc[GRm.index.get_level_values(0) == "4h", "baseMean"] / bm0h

GRm.loc[GRm.index.get_level_values(0) == "8h", "baseMean"] / bm0h

GRm.loc[GRm.index.get_level_values(0) == "12h", "baseMean"] / bm0h





GR["duration"] == "0h"


res0h["baseMean"] / (res1h["baseMean"] + 1)


def mean(l):
    if len(l) > 0:
        return sum(l)/len(l)


finalS = starr[["#'chr'", "'start'", "'end'"]]
finalS["input"] = starr[[f"'input.GSE114063.pool{i}.final.bam'" for i in range(1,13)]].apply(mean, axis=1)
finalS["dex.0h"] = starr[[f"'GR.0h.dex.rep{i}.final.bam'" for i in range(1,6)]].apply(mean, axis=1) / finalS["input"]
finalS["dex.1h"] = starr[[f"'GR.1h.dex.rep{i}.final.bam'" for i in range(1,6)]].apply(mean, axis=1) / finalS["input"]
finalS["dex.4h"] = starr[[f"'GR.4h.dex.rep{i}.final.bam'" for i in range(1,6)]].apply(mean, axis=1) / finalS["input"]
finalS["dex.8h"] = starr[[f"'GR.8h.dex.rep{i}.final.bam'" for i in range(1,6)]].apply(mean, axis=1) / finalS["input"]
finalS["dex.12h"] = starr[[f"'GR.12h.dex.rep{i}.final.bam'" for i in range(1,6)]].apply(mean, axis=1) / finalS["input"]


finalS["dex.0h"]







1.final.bam
