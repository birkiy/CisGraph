


from Functions.Packages import *

# sureCounts = pd.read_table("/groups/lackgrp/ll_members/berkay/enhancerPromoterModel/GSE78709_SuRE-counts-K562_B45_B55_LP170105.txt.gz")



sureCounts = pd.read_table("/groups/lackgrp/ll_members/berkay/enhancerPromoterModel/supplementary_dataset1/SuRE-peaks_K562.45.55_raw_sep_globalLambda.annotated_LP160616.txt.gz")


for i,j  in zip(a,b):
    if i

import pandas as pd
import numpy as np

import seaborn as sns
from matplotlib import pyplot, patches
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import gridspec, colors
from statannot import add_stat_annotation

home="/groups/lackgrp/ll_members/berkay/enhancerPromoterModel"

starrCounts = pd.read_table(home + "/StarrSignal/AR.starr.signal.tsv")

header = ["chr", "start", "end", "name"]
ARBS = pd.read_table(home +"/StarrSignal/ARBS.bed", names=header)

starrCounts = starrCounts.sort_values(["#'chr'", "'start'"]).reset_index(drop=True)
ARBS = ARBS.sort_values(["chr", "start"]).reset_index(drop=True)


starrCounts["name"] = ARBS["name"]


ks =  ["#'chr'","'start'","'end'","'lncap-dht-merged.bam'","'lncap-etoh-merged.bam'","'starrseq-lncap-lib3-cleaned.bam'","name"]
vs = ["chr", "start", "end", "starr.dht", "starr.etoh", "starr.plasmid", "name"]

starrCounts = starrCounts.rename(columns={k:v for k,v in zip(ks, vs)})

starrCounts["normalized.dht"] = np.log(starrCounts["starr.dht"] +1 / starrCounts["starr.plasmid"] +1)
starrCounts["normalized.etoh"] = np.log(starrCounts["starr.etoh"] +1 / starrCounts["starr.plasmid"] +1)


#####################3

home="/groups/lackgrp/ll_members/berkay/enhancerPromoterModel"

header = ["chr", "start", "end", "nameSure", "score", "SuRE.score", "chr.ARBS", "start.ARBS", "end.ARBS", "name"]
con = pd.read_table(home +"/intersect/con.sure.bed", names=header)
con["nodeClass"] = "con"
ind = pd.read_table(home +"/intersect/ind.sure.bed", names=header)
ind["nodeClass"] = "ind"
non = pd.read_table(home +"/intersect/non.sure.bed", names=header)
non["nodeClass"] = "non"
nAR = pd.read_table(home +"/intersect/nAR.sure.bed", names=header)
nAR["nodeClass"] = "nAR"

header = ["chr", "start", "end", "nameSure", "score", "SuRE.score", "chr.ARBS", "start.ARBS", "end.ARBS", "name", "strand"]
pro = pd.read_table(home +"/intersect/pro.sure.bed", names=header)
pro["nodeClass"] = "pro"

DF = pd.concat([con, ind, non, nAR, pro])



mDF = pd.merge(starrCounts, DF, how='inner', on=['name'])



colorPaletteAR = {"con": "#5A5A5A",
                  "ind": "#F9746D",
                  "non": "#ACACAC",
                  "nAR": "#D5D5D5"}

fig = plt.figure(figsize=[16,9])
gs = gridspec.GridSpec(ncols=4, nrows=2)
plt.subplots_adjust(wspace=0.4, hspace=0.4)

params = dict(
    y="SuRE.score"
)

# plt.title(f"DEX treatment for {condition} hours")
for i, nCl in enumerate(["con", "ind", "non", "nAR"]):
    i = i*2
    ax1 = fig.add_subplot(gs[i])
    sns.regplot(**params,x="normalized.etoh", data=mDF[mDF["nodeClass"] == nCl], color=colorPaletteAR[nCl])
    ax1 = fig.add_subplot(gs[i+1])
    sns.regplot(**params,x="normalized.dht", data=mDF[mDF["nodeClass"] == nCl], color=colorPaletteAR[nCl])

fig.savefig(home + "/sureVsStarr.pdf")




g = sns.lmplot(x="normalized.dht", y="SuRE.score", hue="nodeClass", data=mDF, palette=colorPaletteAR)

g.savefig(home + "/dht.sureVsStarr.pdf")


g = sns.lmplot(x="normalized.etoh", y="SuRE.score", hue="nodeClass", data=mDF, palette=colorPaletteAR)

g.savefig(home + "/etoh.sureVsStarr.pdf")










header = ["chr", "start", "end", "nameSure", "score", "SuRE.score", "chr.ARBS", "start.ARBS", "end.ARBS", "name.ARBS", "strand"]
pro = pd.read_table(home +"/pro.sure.bed", names=header)[["chr", "start", "end", "nameSure", "score", "SuRE.score", "chr.ARBS", "start.ARBS", "end.ARBS", "name.ARBS"]]
pro["nodeClass"] = "pro"


DF.sort_values(["chr.ARBS", "start.ARBS"])
starrCounts.sort_values(["#'chr'", "'start'"])



#colorPalette = ["#63b7af", "#abf0e9", "#d4f3ef", "#f5fffd", "#ee8572"]
colorPalette = {"con": "#5A5A5A",
                "ind": "#F9746D",
                "non": "#ACACAC",
                "nAR": "#F5F5F5",
                "pro": "#000000"}


fig = plt.figure(figsize=[8, 5])
gs = gridspec.GridSpec(ncols=2, nrows=1)
plt.subplots_adjust(wspace=0.4)


params = dict(x='nodeClass',
              y='SuRE.score',
              order=["con", "ind", "non", "nAR", "pro"])

boxPairs = [("nAR", "con"),
            ("con", "ind"),
            ("con", "non"),
            ("ind", "non"),
            ("con", "pro")
            ]

ax1 = fig.add_subplot(gs[0])
ax1 = sns.boxplot(**params, data=DF,  palette=colorPalette)

plt.ylabel("SuRE Score", fontsize=14)
plt.xlabel("Node Classes", fontsize=14)
for tick in ax1.xaxis.get_major_ticks():
    tick.label.set_fontsize(13)


add_stat_annotation(ax1, **params, data=DF, box_pairs=boxPairs, test='Mann-Whitney')


fig.savefig(f"{home}/Sure.pdf")














sharpr = pd.read_table("/groups/lackgrp/ll_members/berkay/enhancerPromoterModel/K562_SHARPR-MPRA_scores/basepredictions_K562_ScaleUpDesign1and2_minP.txt", names=["ID"] + list(range(295)))




sud1min1 = pd.read_table("/groups/lackgrp/ll_members/berkay/enhancerPromoterModel/Scaleup_counts_sequences/K562/K562_ScaleUpDesign1_minP_mRNA_Rep1.counts", names=["ID", "min.1"])
sud1min2 = pd.read_table("/groups/lackgrp/ll_members/berkay/enhancerPromoterModel/Scaleup_counts_sequences/K562/K562_ScaleUpDesign1_minP_mRNA_Rep2.counts", names=["ID", "min.2"])

sud1sv1 = pd.read_table("/groups/lackgrp/ll_members/berkay/enhancerPromoterModel/Scaleup_counts_sequences/K562/K562_ScaleUpDesign1_SV40P_mRNA_Rep1.counts", names=["ID", "SV40.1"])
sud1sv2 = pd.read_table("/groups/lackgrp/ll_members/berkay/enhancerPromoterModel/Scaleup_counts_sequences/K562/K562_ScaleUpDesign1_SV40P_mRNA_Rep2.counts", names=["ID", "SV40.2"])

dnamin = pd.read_table("/groups/lackgrp/ll_members/berkay/enhancerPromoterModel/Scaleup_counts_sequences/DNACOUNTS/ScaleUpDesign1_minP_Plasmid.counts", names=["ID", "min.plasmid"])
dnasv = pd.read_table("/groups/lackgrp/ll_members/berkay/enhancerPromoterModel/Scaleup_counts_sequences/DNACOUNTS/ScaleUpDesign1_SV40P_Plasmid.counts", names=["ID", "SV40.plasmid"])




sharprCounts = dnamin
sharprCounts["SV40.plasmid"] = dnasv["SV40.plasmid"]

sharprCounts["min.1"] = sud1min1["min.1"]
sharprCounts["min.2"] = sud1min2["min.2"]

sharprCounts["SV40.1"] = sud1sv1["SV40.1"]
sharprCounts["SV40.2"] = sud1sv2["SV40.2"]

sharprCounts = sharprCounts.drop(0)



new = sharprCounts["ID"].astype(str).str.split("_", expand=True)

sharprCounts["chr"] = new[4]

sharprCounts["center"] = new[5].astype("int64")



sharprCounts = sharprCounts.sort_values("center", ascending=True).reset_index(drop=True)



def find_le(a, x):
    'Find rightmost value less than or equal to x'
    i = bi.bisect_right(a, x)
    if i:
        return i-1
    raise ValueError

def find_ge(a, x):
    'Find leftmost item greater than or equal to x'
    i = bi.bisect_left(a, x)
    if i != len(a):
        return i
    raise ValueError


mapChr = {}
for i in range(BedAll.shape[0]):
    chr = BedAll.loc[i, "chr"]
    center =  BedAll.loc[i, "center"]
    if  not chr in mapChr.keys():
        mapChr[chr] = [int(center)]
    else:
        mapChr[chr].append(int(center))



map2 = {}
for i in range(sureCounts.shape[0]):
    chr = sureCounts.loc[i, "chr"]
    start = sureCounts.loc[i, "start"]
    end = sureCounts.loc[i, "end"]
    s = bi.bisect_left(mapChr[chr], start)
    e = bi.bisect_right(mapChr[chr], end)
    if not chr in map2.keys():
        map2[chr] = mapChr[chr][s:e]
    else:
        map2[chr] += mapChr[chr][s:e]



map2.values()





conBed = pd.read_csv(home + "/ARBSs/regions/cons-arbs.bed", sep="\t", names=["chr", "start", "end", "name"])
conBed["nodeClass"] = "con"
indBed = pd.read_csv(home + "/ARBSs/regions/ind-arbs.bed", sep="\t", names=["chr", "start", "end", "name"])
indBed["nodeClass"] = "ind"
nonBed = pd.read_csv(home + "/ARBSs/regions/Non-Active-ARBS.bed", sep="\t", names=["chr", "start", "end", "name"])
nonBed["nodeClass"] = "non"
nARBed = pd.read_csv(home + "/ARBSs/regions/negativeControl.ARBS.bed", sep="\t", names=["chr", "start", "end", "name"])
nARBed["nodeClass"] = "nAR"
tsPBed = pd.read_csv(home + "/genomeAnnotations/Regions/TSS.hg19.+.bed", sep="\t", names=["chr", "start", "end", "name", "strand"])[["chr", "start", "end", "name"]]
tsPBed["nodeClass"] = "pro"
tsMBed = pd.read_csv(home + "/genomeAnnotations/Regions/TSS.hg19.-.bed", sep="\t", names=["chr", "start", "end", "name", "strand"])[["chr", "start", "end", "name"]]
tsMBed["nodeClass"] = "pro"



BedAll = pd.concat([conBed, indBed, nonBed, nARBed, tsPBed, tsMBed])


BedAll = BedAll.reset_index(drop=True)

pickle.dump(BedAll, open(f"{dataRoot}/BedAll.p", "wb" ))

BedAll.write_table()



BedAll["center"] = ((BedAll["end"] - BedAll["start"])/2).map(int) + BedAll["start"]

BedAll = BedAll.sort_values("center").reset_index(drop=True)
