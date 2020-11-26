




import numpy as np
import pandas as pd


import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import gridspec, colors
# from pygraphviz import *
from statannot import add_stat_annotation


header = [
  "chr", "start", "end",
  "GR.0h.dex.rep1", "GR.0h.dex.rep2", "GR.0h.dex.rep3", "GR.0h.dex.rep4", "GR.0h.dex.rep5",
  "GR.12h.dex.rep1", "GR.12h.dex.rep2", "GR.12h.dex.rep3", "GR.12h.dex.rep4", "GR.12h.dex.rep5",
  "GR.1h.dex.rep1", "GR.1h.dex.rep2", "GR.1h.dex.rep3", "GR.1h.dex.rep4", "GR.1h.dex.rep5",
  "GR.4h.dex.rep1", "GR.4h.dex.rep2", "GR.4h.dex.rep3", "GR.4h.dex.rep4", "GR.4h.dex.rep5",
  "GR.8h.dex.rep1", "GR.8h.dex.rep2", "GR.8h.dex.rep3", "GR.8h.dex.rep4", "GR.8h.dex.rep5",
  "input.GSE114063.pool10", "input.GSE114063.pool11","input.GSE114063.pool12", "input.GSE114063.pool1", "input.GSE114063.pool2", "input.GSE114063.pool3", "input.GSE114063.pool4","input.GSE114063.pool5","input.GSE114063.pool6","input.GSE114063.pool7","input.GSE114063.pool8","input.GSE114063.pool9"

]

countTable = pd.read_table("/groups/lackgrp/ll_members/berkay/StarrPipe/results/coverage/countTableRaw.final.txt", names=header, skiprows=1)


home = "/groups/lackgrp/ll_members/berkay/STARRbegin/peaks"

bedGR = pd.read_table(f"{home}/common.GR.peaks.bed", names=["chr", "start", "end", "name", "score"])
bedNGR = pd.read_table(f"{home}/negativeControlGR.final.bed", names=["chr", "start", "end", "name"])


def posName(bed):
    return bed["chr"] +":"+ bed["start"].map(str) +"-"+ bed["end"].map(str)


bedNGR["posName"] = posName(bedNGR)

bedGR["posName"] = posName(bedGR)

countTable["posName"] = posName(countTable)


countTable.loc[countTable["posName"].isin(bedNGR["posName"]), "nodeClass"] = "nGR"
countTable.loc[countTable["posName"].isin(bedGR["posName"]), "nodeClass"] = "GR"



h0 = countTable[["GR.0h.dex.rep1", "GR.0h.dex.rep2", "GR.0h.dex.rep3", "GR.0h.dex.rep4", "GR.0h.dex.rep5"]].values
h8 = countTable[["GR.8h.dex.rep1", "GR.8h.dex.rep2", "GR.8h.dex.rep3", "GR.8h.dex.rep4", "GR.8h.dex.rep5"]].values

inp = countTable[["input.GSE114063.pool10", "input.GSE114063.pool11","input.GSE114063.pool12", "input.GSE114063.pool1", "input.GSE114063.pool2", "input.GSE114063.pool3", "input.GSE114063.pool4","input.GSE114063.pool5", "input.GSE114063.pool6","input.GSE114063.pool7","input.GSE114063.pool8","input.GSE114063.pool9"]].values



def batchNormalization(h0, e=0.001):
  μ = h0.mean()
  σ = h0.var()
  ϵ = e
  return (h0 - μ) / (σ - ϵ)**0.5




h0 = batchNormalization(h0)
h8 = batchNormalization(h8)
inp = batchNormalization(inp)


df = pd.DataFrame(np.concatenate((h0,h8,inp), axis=1), columns=["GR.0h.dex.rep1", "GR.0h.dex.rep2", "GR.0h.dex.rep3", "GR.0h.dex.rep4", "GR.0h.dex.rep5", "GR.8h.dex.rep1", "GR.8h.dex.rep2", "GR.8h.dex.rep3", "GR.8h.dex.rep4", "GR.8h.dex.rep5", "input.GSE114063.pool10", "input.GSE114063.pool11","input.GSE114063.pool12", "input.GSE114063.pool1", "input.GSE114063.pool2", "input.GSE114063.pool3", "input.GSE114063.pool4","input.GSE114063.pool5","input.GSE114063.pool6","input.GSE114063.pool7","input.GSE114063.pool8","input.GSE114063.pool9"])

normCounts = pd.concat([countTable[["chr", "start", "end", "posName", "nodeClass"]], df], axis=1)




sumCount = countTable[["chr", "start", "end", "posName", "nodeClass"]]
sumCount["0h"] = normCounts[["GR.0h.dex.rep1", "GR.0h.dex.rep2", "GR.0h.dex.rep3", "GR.0h.dex.rep4", "GR.0h.dex.rep5"]].sum(axis=1)
sumCount["8h"] = normCounts[["GR.8h.dex.rep1", "GR.8h.dex.rep2", "GR.8h.dex.rep3", "GR.8h.dex.rep4", "GR.8h.dex.rep5"]].sum(axis=1)
sumCount["input"] = normCounts[["input.GSE114063.pool10", "input.GSE114063.pool11","input.GSE114063.pool12", "input.GSE114063.pool1", "input.GSE114063.pool2", "input.GSE114063.pool3", "input.GSE114063.pool4","input.GSE114063.pool5","input.GSE114063.pool6","input.GSE114063.pool7","input.GSE114063.pool8","input.GSE114063.pool9"]].sum(axis=1)


sumCount["plasmid.0h"] = np.log2((sumCount["0h"] +1) / (sumCount["input"] +1))
sumCount["plasmid.8h"] = np.log2((sumCount["8h"] +1) / (sumCount["input"] +1))

n = sumCount[(sumCount["plasmid.8h"] > 1) & (sumCount["plasmid.0h"] > 1)]

pd.concat([countTable[["chr", "start", "end", "posName", "nodeClass"]], n ], axis=1, ignore_axis=True)















#####################

sumCount = pd.DataFrame()
sumCount["0h"] = countTable[["GR.0h.dex.rep1", "GR.0h.dex.rep2", "GR.0h.dex.rep3", "GR.0h.dex.rep4", "GR.0h.dex.rep5"]].sum(axis=1)
sumCount["8h"] = countTable[["GR.8h.dex.rep1", "GR.8h.dex.rep2", "GR.8h.dex.rep3", "GR.8h.dex.rep4", "GR.8h.dex.rep5"]].sum(axis=1)
sumCount["input"] = countTable[["input.GSE114063.pool10", "input.GSE114063.pool11","input.GSE114063.pool12", "input.GSE114063.pool1", "input.GSE114063.pool2", "input.GSE114063.pool3", "input.GSE114063.pool4","input.GSE114063.pool5","input.GSE114063.pool6","input.GSE114063.pool7","input.GSE114063.pool8","input.GSE114063.pool9"]].sum(axis=1)


def quantileNormalization(df):
  rankMean = df.stack().groupby(df.rank(method='first').stack().astype(int)).mean()
  df = df.rank(method='min').stack().astype(int).map(rankMean).unstack()
  return df


df = quantileNormalization(sumCount[["0h", "8h", "input"]])

dfm = df.melt()


g = sns.displot(dfm, x="value", hue="variable")

g.savefig("/groups/lackgrp/ll_members/berkay/STARRbegin/0h.8h.lib.qn.pdf")



sumCount["0h"] = sumCount["0h"]
sumCount["8h"] = sumCount["8h"]
sumCount["input"] = sumCount["input"] / (sumCount["input"].sum() ** (1/2))


df["plasmid.0h"] = np.log2((df["0h"] +1) / (df["input"] +1))
df["plasmid.8h"] = np.log2((df["8h"] +1) / (df["input"] +1))





sumCount[(sumCount["plasmid.8h"] > 1) & (sumCount["plasmid.0h"] > 1)]









sumCount = pd.DataFrame()
sumCount["0h"] = countTable[["GR.0h.dex.rep1", "GR.0h.dex.rep2", "GR.0h.dex.rep3", "GR.0h.dex.rep4", "GR.0h.dex.rep5"]].sum(axis=1)
sumCount["8h"] = countTable[["GR.8h.dex.rep1", "GR.8h.dex.rep2", "GR.8h.dex.rep3", "GR.8h.dex.rep4", "GR.8h.dex.rep5"]].sum(axis=1)
sumCount["input"] = countTable[["input.GSE114063.pool10", "input.GSE114063.pool11","input.GSE114063.pool12", "input.GSE114063.pool1", "input.GSE114063.pool2", "input.GSE114063.pool3", "input.GSE114063.pool4","input.GSE114063.pool5","input.GSE114063.pool6","input.GSE114063.pool7","input.GSE114063.pool8","input.GSE114063.pool9"]].mean(axis=1)

sumCount["plasmid.0h"] = np.log2((sumCount["0h"] +1) / (sumCount["input"] +1))
sumCount["plasmid.8h"] = np.log2((sumCount["8h"] +1) / (sumCount["input"] +1))


sumCount["nodeClass"] = countTable["nodeClass"]



sumCount["posName"] = countTable["posName"]



sumCount["plasmid.8h"] / sumCount["plasmid.0h"]




sumCount["0h"] = countTable[["GR.0h.dex.rep1", "GR.0h.dex.rep2", "GR.0h.dex.rep3", "GR.0h.dex.rep4", "GR.0h.dex.rep5"]].sum(axis=1)
sumCount["0h"] = countTable[["GR.0h.dex.rep1", "GR.0h.dex.rep2", "GR.0h.dex.rep3", "GR.0h.dex.rep4", "GR.0h.dex.rep5"]].sum(axis=1)


##################################################################



def readDESeq2(file, duration):
    DF = pd.read_table(file)
    new = DF.reset_index()["index"].str.split("-", expand=True)
    DF = DF.reset_index(drop=True)
    DF["chr"] = new[0].str.strip()
    DF["start"] = new[1].str.strip().astype("int64")
    DF["end"] = new[2].str.strip().astype("int64")
    DF["duration"] = duration
    return DF


lib0 = readDESeq2("/groups/lackgrp/ll_members/berkay/StarrPipe/results/coverage/GR.0h.lib.tsv", "0h")
lib1 = readDESeq2("/groups/lackgrp/ll_members/berkay/StarrPipe/results/coverage/GR.1h.lib.tsv", "1h")
lib4 = readDESeq2("/groups/lackgrp/ll_members/berkay/StarrPipe/results/coverage/GR.4h.lib.tsv", "4h")
lib8 = readDESeq2("/groups/lackgrp/ll_members/berkay/StarrPipe/results/coverage/GR.8h.lib.tsv", "8h")
lib12 = readDESeq2("/groups/lackgrp/ll_members/berkay/StarrPipe/results/coverage/GR.12h.lib.tsv", "12h")

lib = pd.concat([lib0, lib1, lib4, lib8, lib12]).reset_index(drop=True)
lib["comparison"] = "lib"


home = "/groups/lackgrp/ll_members/berkay/STARRbegin/peaks"

bedGR = pd.read_table(f"{home}/common.GR.peaks.bed", names=["chr", "start", "end", "name", "score"])
bedNGR = pd.read_table(f"{home}/negativeControlGR.final.bed", names=["chr", "start", "end", "name"])

def posName(bed):
    return bed["chr"] +":"+ bed["start"].map(str) +"-"+ bed["end"].map(str)

bedNGR["posName"] = posName(bedNGR)
bedGR["posName"] = posName(bedGR)
lib["posName"] = posName(lib)

lib.loc[lib["posName"].isin(bedNGR["posName"]), "nodeClass"] = "nGR"
lib.loc[lib["posName"].isin(bedGR["posName"]), "nodeClass"] = "GR"


z = [
    (lfc, dur)
    for lfc in [0.5, 0.6, 0.7, 0.8, 0.9, 1]
    for dur in ["1h", "4h", "8h", "12h"]
]

for lfc, dur in z:
    ind = pd.read_table(f"/groups/lackgrp/ll_members/berkay/StarrPipe/regions/ind.{lfc}.GR.{dur}.bed", names=["chr", "start", "end", "posName", "nodeClass"])
    ind["posName"] = posName(ind)
    libSub = lib[~(lib["posName"].isin(ind["posName"]))]
    con = libSub.loc[
        (libSub["log2FoldChange"] >= lfc) &
        (libSub["nodeClass"] == "GR") &
        (libSub["padj"] < 0.05) &
        (libSub["duration"] == dur),
        ["chr", "start", "end", "posName", "nodeClass"]
    ]
    print(f"CON => lfc:{lfc}\tdur:{dur}\tn:{len(con)}")
    con["nodeClass"] = "con"
    non = libSub.loc[
        ~(libSub["posName"].isin(list(ind["posName"]) + list(con["posName"]))) &
        (libSub["nodeClass"] == "GR") &
        (libSub["duration"] == dur),
        ["chr", "start", "end","posName", "nodeClass"]
    ]
    print(f"NON => lfc:{lfc}\tdur:{dur}\tn:{len(non)}")
    non["nodeClass"] = "non"
    con.to_csv(f"/groups/lackgrp/ll_members/berkay/StarrPipe/regions/con.{lfc}.GR.{dur}.bed",sep="\t", index=False, header=False)
    non.to_csv(f"/groups/lackgrp/ll_members/berkay/StarrPipe/regions/non.{lfc}.GR.{dur}.bed",sep="\t", index=False, header=False)







libSub = libSub[libSub["log2FoldChange"] > 0]

con2 = pd.read_table("/groups/lackgrp/ll_members/berkay/STARRbegin/peaks/con.GR.bed", names=["chr", "start", "end", "posName", "nodeClass"])



["chr", "start", "end","posName", "nodeClass", "padj", "log2FoldChange"]
########################3

libt = readDESeq2("/groups/lackgrp/ll_members/berkay/STARRbegin/results/coverage/GR.starr.lib.tsv", "starr")

libt["comparison"] = "lib"
libt = libt[libt["padj"] < 0.05]
libt["posName"] = posName(libt)

libt.loc[libt["posName"].isin(bedNGR["posName"]), "nodeClass"] = "nGR"
libt.loc[libt["posName"].isin(bedGR["posName"]), "nodeClass"] = "GR"

cont = libt.loc[(libt["log2FoldChange"] >= 1) & (libt["nodeClass"] == "GR"),:].drop_duplicates()
 ["chr", "start", "end","posName", "nodeClass"]
con = pd.concat([con,cont]).drop_duplicates().sort_values(["chr", "start"])
con["nodeClass"] = "con"

non = libSub.loc[~(libSub["posName"].isin(list(ind["posName"]) + list(con["posName"]))), ["chr", "start", "end","posName", "nodeClass"]].drop_duplicates().sort_values(["chr", "start"])
non["nodeClass"] = "non"



con.to_csv("/groups/lackgrp/ll_members/berkay/STARRbegin/peaks/con.GR.bed",sep="\t", index=False, header=False)

non.to_csv("/groups/lackgrp/ll_members/berkay/STARRbegin/peaks/non.GR.bed",sep="\t", index=False, header=False)
