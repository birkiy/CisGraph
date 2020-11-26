





import pandas as pd

import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import ticker as mticker
from statannot import add_stat_annotation

from sklearn.cluster import KMeans



metadata = pd.read_table("/groups/lackgrp/ll_members/berkay/DHS/rawData/DHS_Index_and_Vocabulary_metadata.tsv")
samples = metadata[metadata["System"] == "Genitourinary"]["DCC File ID"].values

indexARBS = pd.read_table("/groups/lackgrp/ll_members/berkay/DHS/creOverARBS/ARBS.index.txt", names=["index"])

binarydata = np.loadtxt("/groups/lackgrp/ll_members/berkay/DHS/creOverARBS/indexed.binary.txt", delimiter="\t")



systems = metadata.System.unique()[:-1]
systemDict = {}
for system in systems:
    idxSamples = metadata[metadata["System"] == system].index
    subMat = np.array([binarydata[:,i] for i in idxSamples])
    systemDict[system] = (subMat.sum(axis=0) / len(idxSamples)) * 100

systemDF = pd.DataFrame(systemDict, index=indexARBS["index"])

idxSamples = metadata[metadata["System"] == "Genitourinary"].index
subMat = np.array([binarydata[:,i] for i in idxSamples])

binDF = pd.DataFrame(data=np.transpose(subMat), columns=samples, index=indexARBS["index"])




header = ["chr", "start", "end", "identifier", "mean_signal", "numsamples", "index.ARBS"]
con = pd.read_table("/groups/lackgrp/ll_members/berkay/DHS/creOverARBS/con.cre.bed", names=header)
con["nodeClass"] = "con"
ind = pd.read_table("/groups/lackgrp/ll_members/berkay/DHS/creOverARBS/ind.cre.bed", names=header)
ind["nodeClass"] = "ind"
non = pd.read_table("/groups/lackgrp/ll_members/berkay/DHS/creOverARBS/non.cre.bed", names=header)
non["nodeClass"] = "non"
nAR = pd.read_table("/groups/lackgrp/ll_members/berkay/DHS/creOverARBS/nAR.cre.bed", names=header)
nAR["nodeClass"] = "nAR"




cre = pd.concat([con, ind, non, nAR]).reset_index(drop=True)
cre["nlog"] = np.log10(cre["numsamples"])

systemDict = systemDF.to_dict(orient="index")
binDict = binDF.to_dict(orient="index")
creDict = cre.to_dict(orient="index")



concatDict = {}
for i in range(cre.shape[0]):
    idx = cre.loc[i, "index.ARBS"]
    tmp = dict(**creDict[i], **binDict[idx], **systemDict[idx])
    concatDict[i] = tmp

creBin = pd.DataFrame.from_dict(concatDict, orient="index")

creBin = creBin.sort_values(["nodeClass", "Genitourinary"], ascending=False)


####################################################



fig = plt.figure(figsize=[5,6])

colorPaletteAR = {"con": "#5A5A5A",
                  "ind": "#F9746D",
                  "non": "#ACACAC",
                  "nAR": "#D5D5D5"}

boxPairs = [("con", "ind"),
            ("con", "non"),
            ("ind", "non")
            ]
ax = sns.boxplot(data=cre, y="nlog", x="nodeClass", palette = colorPaletteAR ,linewidth=2)
#ax = sns.boxplot(data=cre, y="nlog", x="nodeClass", color = "black", width = .15, zorder = 10, showcaps = False, boxprops = {'facecolor':'#FFFFFF', "zorder":10}, showfliers=False, whiskerprops = {'linewidth':2, "zorder":10}, saturation = 1)
plt.xlabel("NodeClasses", fontsize=14)
plt.ylabel("# samples CRE found in", fontsize=14)
ax.yaxis.set_ticks([0,1,2,3])
ax.yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
ax.yaxis.set_ticks([np.log10(x) for p in range(0,3) for x in np.linspace(10**p, 10**(p+1), 10)], minor=True)
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(13)

for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(13)

add_stat_annotation(ax, data=cre, y="nlog", x="nodeClass", box_pairs=boxPairs, test='Mann-Whitney', loc='inside' ,line_height=0, text_offset=1,line_offset=0.005)

fig.tight_layout()
fig.savefig("/groups/lackgrp/ll_members/berkay/DHS/creOverARBS/creARBS.error.pdf")



###############################################


fig = plt.figure(figsize=[5,5])

colorPaletteAR = {"con": "#5A5A5A",
                  "ind": "#F9746D",
                  "non": "#ACACAC",
                  "nAR": "#D5D5D5"}

boxPairs = [("con", "ind"),
            ("con", "non"),
            ("ind", "non"),
            ("nAR", "ind"),
            ("con", "nAR"),
            ("nAR", "non")
            ]
ax = sns.stripplot(data=cre, y="numsamples", x="nodeClass", palette = colorPaletteAR)
plt.xlabel("NodeClasses", fontsize=14)
plt.ylabel("# samples CRE found in", fontsize=14)
add_stat_annotation(ax, data=cre, y="numsamples", x="nodeClass", box_pairs=boxPairs, test='Mann-Whitney')

fig.tight_layout()
fig.savefig("/groups/lackgrp/ll_members/berkay/DHS/creOverARBS/creARBS.strip.pdf")

###############################################




prostateTissues = [
    'ENCFF503PAE', 'ENCFF782HVG',
    'ENCFF939VCZ', 'ENCFF329UDK', 'ENCFF015ICH', 'ENCFF658UUI',
    'ENCFF335MSP', 'ENCFF537ICV', 'ENCFF081ETZ', 'ENCFF875DXE',
    'ENCFF898KEE', 'ENCFF309DYE', 'ENCFF176OMF', 'ENCFF625QKW'
    ]
labels = ['Hela', 'HeLaS3', 'PrEC', 'PrEC', 'LNCap', 'LNCap', 'fTestes',
       'fTestes', 'Ovary', 'fOvary', 'HeLaS3', 'HeLaS3', 'PC3', 'PC3']

creBin = creBin.rename(columns=dict(zip(prostateTissues, labels)))


specificProstate = ['PrEC', 'PrEC', 'LNCap', 'LNCap', 'PC3', 'PC3']
creBin[specificProstate]


row_colors = creBin["nodeClass"].map(colorPaletteAR)
cm = sns.light_palette("#69d", as_cmap=True)
g = sns.clustermap(creBin[specificProstate], row_colors=row_colors, cmap=cm, row_cluster=False)

g.savefig("/groups/lackgrp/ll_members/berkay/DHS/creOverARBS/ProstateTissues.pdf")


###############################################
