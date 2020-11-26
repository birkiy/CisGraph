

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

countTable = pd.read_table("/groups/lackgrp/ll_members/berkay/STARRbegin/results/coverage/countTableRaw.txt", names=header, skiprows=1)


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


lib0 = readDESeq2("/groups/lackgrp/ll_members/berkay/STARRbegin/results/coverage/GR.0h.lib.tsv", "0h")
lib8 = readDESeq2("/groups/lackgrp/ll_members/berkay/STARRbegin/results/coverage/GR.8h.lib.tsv", "8h")

lib = pd.concat([lib0, lib8]).reset_index(drop=True)
lib["comparison"] = "lib"
lib = lib[lib["padj"] < 0.05]


ind["posName"] = posName(ind)
lib["posName"] = posName(lib)

lib.loc[lib["posName"].isin(bedNGR["posName"]), "nodeClass"] = "nGR"
lib.loc[lib["posName"].isin(bedGR["posName"]), "nodeClass"] = "GR"

libSub = lib[~(lib["posName"].isin(ind["posName"]))]

libSub = libSub[libSub["log2FoldChange"] > 0]

con = libSub.loc[(libSub["log2FoldChange"] >= 1) & (libSub["nodeClass"] == "GR"), ["chr", "start", "end","posName", "nodeClass"]].drop_duplicates()

########################3

libt = readDESeq2("/groups/lackgrp/ll_members/berkay/STARRbegin/results/coverage/GR.starr.lib.tsv", "starr")

libt["comparison"] = "lib"
libt = libt[libt["padj"] < 0.05]
libt["posName"] = posName(libt)

libt.loc[libt["posName"].isin(bedNGR["posName"]), "nodeClass"] = "nGR"
libt.loc[libt["posName"].isin(bedGR["posName"]), "nodeClass"] = "GR"

cont = libt.loc[(libt["log2FoldChange"] >= 1) & (libt["nodeClass"] == "GR"), ["chr", "start", "end","posName", "nodeClass"]].drop_duplicates()

con = pd.concat([con,cont]).drop_duplicates().sort_values(["chr", "start"])
con["nodeClass"] = "con"

non = libSub.loc[~(libSub["posName"].isin(list(ind["posName"]) + list(con["posName"]))), ["chr", "start", "end","posName", "nodeClass"]].drop_duplicates().sort_values(["chr", "start"])
non["nodeClass"] = "non"



con.to_csv("/groups/lackgrp/ll_members/berkay/STARRbegin/peaks/con.GR.bed",sep="\t", index=False, header=False)

non.to_csv("/groups/lackgrp/ll_members/berkay/STARRbegin/peaks/non.GR.bed",sep="\t", index=False, header=False)
