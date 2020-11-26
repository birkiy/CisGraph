

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from statannot import add_stat_annotation
from matplotlib import gridspec, colors
from matplotlib import ticker as mticker
import numpy as np

home = "/groups/lackgrp/ll_members/berkay/enhancerPromoterModel/SuREsnpCount/ARBS"

header = ["chr", "start", "end", "strand", "lib", "K562.1", "K562.2", "K562.3", "HepG2.1", "HepG2.2", "chr.ARBS", "start.ARBS", "end.ARBS", "name.ARBS"]


SuRE421con = pd.read_table(f"{home}/GSE128325_SuRE42_1.count.con.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE421con["individual"] = "SuRE42"
SuRE421con["rep"] = "1"
SuRE421con["nodeClass"] = "con"
SuRE421ind = pd.read_table(f"{home}/GSE128325_SuRE42_1.count.ind.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE421ind["individual"] = "SuRE42"
SuRE421ind["rep"] = "1"
SuRE421ind["nodeClass"] = "ind"
SuRE421nAR = pd.read_table(f"{home}/GSE128325_SuRE42_1.count.nAR.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE421nAR["individual"] = "SuRE42"
SuRE421nAR["rep"] = "1"
SuRE421nAR["nodeClass"] = "nAR"
SuRE421non = pd.read_table(f"{home}/GSE128325_SuRE42_1.count.non.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE421non["individual"] = "SuRE42"
SuRE421non["rep"] = "1"
SuRE421non["nodeClass"] = "non"

SuRE422con = pd.read_table(f"{home}/GSE128325_SuRE42_2.count.con.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE422con["individual"] = "SuRE42"
SuRE422con["rep"] = "2"
SuRE422con["nodeClass"] = "con"
SuRE422ind = pd.read_table(f"{home}/GSE128325_SuRE42_2.count.ind.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE422ind["individual"] = "SuRE42"
SuRE422ind["rep"] = "2"
SuRE422ind["nodeClass"] = "ind"
SuRE422nAR = pd.read_table(f"{home}/GSE128325_SuRE42_2.count.nAR.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE422nAR["individual"] = "SuRE42"
SuRE422nAR["rep"] = "2"
SuRE422nAR["nodeClass"] = "nAR"
SuRE422non = pd.read_table(f"{home}/GSE128325_SuRE42_2.count.non.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE422non["individual"] = "SuRE42"
SuRE422non["rep"] = "2"
SuRE422non["nodeClass"] = "non"

SuRE431con = pd.read_table(f"{home}/GSE128325_SuRE43_1.count.con.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE431con["individual"] = "SuRE43"
SuRE431con["rep"] = "1"
SuRE431con["nodeClass"] = "con"
SuRE431ind = pd.read_table(f"{home}/GSE128325_SuRE43_1.count.ind.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE431ind["individual"] = "SuRE43"
SuRE431ind["rep"] = "1"
SuRE431ind["nodeClass"] = "ind"
SuRE431nAR = pd.read_table(f"{home}/GSE128325_SuRE43_1.count.nAR.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE431nAR["individual"] = "SuRE43"
SuRE431nAR["rep"] = "1"
SuRE431nAR["nodeClass"] = "nAR"
SuRE431non = pd.read_table(f"{home}/GSE128325_SuRE43_1.count.non.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE431non["individual"] = "SuRE43"
SuRE431non["rep"] = "1"
SuRE431non["nodeClass"] = "non"

SuRE432con = pd.read_table(f"{home}/GSE128325_SuRE43_2.count.con.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE432con["individual"] = "SuRE43"
SuRE432con["rep"] = "2"
SuRE432con["nodeClass"] = "con"
SuRE432ind = pd.read_table(f"{home}/GSE128325_SuRE43_2.count.ind.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE432ind["individual"] = "SuRE43"
SuRE432ind["rep"] = "2"
SuRE432ind["nodeClass"] = "ind"
SuRE432nAR = pd.read_table(f"{home}/GSE128325_SuRE43_2.count.nAR.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE432nAR["individual"] = "SuRE43"
SuRE432nAR["rep"] = "2"
SuRE432nAR["nodeClass"] = "nAR"
SuRE432non = pd.read_table(f"{home}/GSE128325_SuRE43_2.count.non.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE432non["individual"] = "SuRE43"
SuRE432non["rep"] = "2"
SuRE432non["nodeClass"] = "non"

SuRE441con = pd.read_table(f"{home}/GSE128325_SuRE44_1.count.con.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE441con["individual"] = "SuRE44"
SuRE441con["rep"] = "1"
SuRE441con["nodeClass"] = "con"
SuRE441ind = pd.read_table(f"{home}/GSE128325_SuRE44_1.count.ind.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE441ind["individual"] = "SuRE44"
SuRE441ind["rep"] = "1"
SuRE441ind["nodeClass"] = "ind"
SuRE441nAR = pd.read_table(f"{home}/GSE128325_SuRE44_1.count.nAR.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE441nAR["individual"] = "SuRE44"
SuRE441nAR["rep"] = "1"
SuRE441nAR["nodeClass"] = "nAR"
SuRE441non = pd.read_table(f"{home}/GSE128325_SuRE44_1.count.non.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE441non["individual"] = "SuRE44"
SuRE441non["rep"] = "1"
SuRE441non["nodeClass"] = "non"

SuRE442con = pd.read_table(f"{home}/GSE128325_SuRE44_2.count.con.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE442con["individual"] = "SuRE44"
SuRE442con["rep"] = "2"
SuRE442con["nodeClass"] = "con"
SuRE442ind = pd.read_table(f"{home}/GSE128325_SuRE44_2.count.ind.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE442ind["individual"] = "SuRE44"
SuRE442ind["rep"] = "2"
SuRE442ind["nodeClass"] = "ind"
SuRE442nAR = pd.read_table(f"{home}/GSE128325_SuRE44_2.count.nAR.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE442nAR["individual"] = "SuRE44"
SuRE442nAR["rep"] = "2"
SuRE442nAR["nodeClass"] = "nAR"
SuRE442non = pd.read_table(f"{home}/GSE128325_SuRE44_2.count.non.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE442non["individual"] = "SuRE44"
SuRE442non["rep"] = "2"
SuRE442non["nodeClass"] = "non"

SuRE451con = pd.read_table(f"{home}/GSE128325_SuRE45_1.count.con.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE451con["individual"] = "SuRE45"
SuRE451con["rep"] = "1"
SuRE451con["nodeClass"] = "con"
SuRE451ind = pd.read_table(f"{home}/GSE128325_SuRE45_1.count.ind.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE451ind["individual"] = "SuRE45"
SuRE451ind["rep"] = "1"
SuRE451ind["nodeClass"] = "ind"
SuRE451nAR = pd.read_table(f"{home}/GSE128325_SuRE45_1.count.nAR.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE451nAR["individual"] = "SuRE45"
SuRE451nAR["rep"] = "1"
SuRE451nAR["nodeClass"] = "nAR"
SuRE451non = pd.read_table(f"{home}/GSE128325_SuRE45_1.count.non.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE451non["individual"] = "SuRE45"
SuRE451non["rep"] = "1"
SuRE451non["nodeClass"] = "non"

SuRE452con = pd.read_table(f"{home}/GSE128325_SuRE45_2.count.con.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE452con["individual"] = "SuRE45"
SuRE452con["rep"] = "2"
SuRE452con["nodeClass"] = "con"
SuRE452ind = pd.read_table(f"{home}/GSE128325_SuRE45_2.count.ind.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE452ind["individual"] = "SuRE45"
SuRE452ind["rep"] = "2"
SuRE452ind["nodeClass"] = "ind"
SuRE452nAR = pd.read_table(f"{home}/GSE128325_SuRE45_2.count.nAR.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE452nAR["individual"] = "SuRE45"
SuRE452nAR["rep"] = "2"
SuRE452nAR["nodeClass"] = "nAR"
SuRE452non = pd.read_table(f"{home}/GSE128325_SuRE45_2.count.non.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE452non["individual"] = "SuRE45"
SuRE452non["rep"] = "2"
SuRE452non["nodeClass"] = "non"







############


header = ["chr", "start", "end", "strand", "lib", "K562.1", "K562.2", "K562.3", "HepG2.1", "HepG2.2", "chr.ARBS", "start.ARBS", "end.ARBS", "name.ARBS", "strand.ARBS"]

SuRE421tss = pd.read_table(f"{home}/GSE128325_SuRE42_1.count.tss.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE421tss["individual"] = "SuRE42"
SuRE421tss["rep"] = "1"
SuRE421tss["nodeClass"] = "pro"

SuRE422tss = pd.read_table(f"{home}/GSE128325_SuRE42_2.count.tss.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE422tss["individual"] = "SuRE42"
SuRE422tss["rep"] = "2"
SuRE422tss["nodeClass"] = "pro"

SuRE431tss = pd.read_table(f"{home}/GSE128325_SuRE43_1.count.tss.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE431tss["individual"] = "SuRE43"
SuRE431tss["rep"] = "1"
SuRE431tss["nodeClass"] = "pro"

SuRE432tss = pd.read_table(f"{home}/GSE128325_SuRE43_2.count.tss.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE432tss["individual"] = "SuRE43"
SuRE432tss["rep"] = "2"
SuRE432tss["nodeClass"] = "pro"

SuRE441tss = pd.read_table(f"{home}/GSE128325_SuRE44_1.count.tss.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE441tss["individual"] = "SuRE44"
SuRE441tss["rep"] = "1"
SuRE441tss["nodeClass"] = "pro"

SuRE442tss = pd.read_table(f"{home}/GSE128325_SuRE44_2.count.tss.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE442tss["individual"] = "SuRE44"
SuRE442tss["rep"] = "2"
SuRE442tss["nodeClass"] = "pro"

SuRE451tss = pd.read_table(f"{home}/GSE128325_SuRE45_1.count.tss.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE451tss["individual"] = "SuRE45"
SuRE451tss["rep"] = "1"
SuRE451tss["nodeClass"] = "pro"

SuRE452tss = pd.read_table(f"{home}/GSE128325_SuRE45_2.count.tss.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
SuRE452tss["individual"] = "SuRE45"
SuRE452tss["rep"] = "2"
SuRE452tss["nodeClass"] = "pro"



SuRE42 = pd.concat([SuRE421con, SuRE421ind, SuRE421nAR, SuRE421non, SuRE421tss, SuRE422con, SuRE422ind, SuRE422nAR, SuRE422non, SuRE422tss]).reset_index(drop=True)
# SuRE42["lib.n"] = SuRE42["lib"] / SuRE42["lib"].sum()
# SuRE42["K562.1.n"] = SuRE42["K562.1"] / SuRE42["K562.1"].sum()
# SuRE42["K562.2.n"] = SuRE42["K562.2"] / SuRE42["K562.2"].sum()
# SuRE42["K562.3.n"] = SuRE42["K562.3"] / SuRE42["K562.3"].sum()
# SuRE42["HepG2.1.n"] = SuRE42["HepG2.1"] / SuRE42["HepG2.1"].sum()
# SuRE42["HepG2.2.n"] = SuRE42["HepG2.2"] / SuRE42["HepG2.2"].sum()

SuRE43 = pd.concat([SuRE431con, SuRE431ind, SuRE431nAR, SuRE431non, SuRE431tss, SuRE432con, SuRE432ind, SuRE432nAR, SuRE432non, SuRE432tss]).reset_index(drop=True)
# SuRE43["lib.n"] = SuRE43["lib"] / SuRE43["lib"].sum()
# SuRE43["K562.1.n"] = SuRE43["K562.1"] / SuRE43["K562.1"].sum()
# SuRE43["K562.2.n"] = SuRE43["K562.2"] / SuRE43["K562.2"].sum()
# SuRE43["K562.3.n"] = SuRE43["K562.3"] / SuRE43["K562.3"].sum()
# SuRE43["HepG2.1.n"] = SuRE43["HepG2.1"] / SuRE43["HepG2.1"].sum()
# SuRE43["HepG2.2.n"] = SuRE43["HepG2.2"] / SuRE43["HepG2.2"].sum()

SuRE44 = pd.concat([SuRE441con, SuRE441ind, SuRE441nAR, SuRE441non, SuRE441tss, SuRE442con, SuRE442ind, SuRE442nAR, SuRE442non, SuRE442tss]).reset_index(drop=True)
# SuRE44["lib.n"] = SuRE44["lib"] / SuRE44["lib"].sum()
# SuRE44["K562.1.n"] = SuRE44["K562.1"] / SuRE44["K562.1"].sum()
# SuRE44["K562.2.n"] = SuRE44["K562.2"] / SuRE44["K562.2"].sum()
# SuRE44["K562.3.n"] = SuRE44["K562.3"] / SuRE44["K562.3"].sum()
# SuRE44["HepG2.1.n"] = SuRE44["HepG2.1"] / SuRE44["HepG2.1"].sum()
# SuRE44["HepG2.2.n"] = SuRE44["HepG2.2"] / SuRE44["HepG2.2"].sum()

SuRE45 = pd.concat([SuRE451con, SuRE451ind, SuRE451nAR, SuRE451non, SuRE451tss, SuRE452con, SuRE452ind, SuRE452nAR, SuRE452non, SuRE452tss]).reset_index(drop=True)
# SuRE45["lib.n"] = SuRE45["lib"] / SuRE45["lib"].sum()
# SuRE45["K562.1.n"] = SuRE45["K562.1"] / SuRE45["K562.1"].sum()
# SuRE45["K562.2.n"] = SuRE45["K562.2"] / SuRE45["K562.2"].sum()
# SuRE45["K562.3.n"] = SuRE45["K562.3"] / SuRE45["K562.3"].sum()
# SuRE45["HepG2.1.n"] = SuRE45["HepG2.1"] / SuRE45["HepG2.1"].sum()
# SuRE45["HepG2.2.n"] = SuRE45["HepG2.2"] / SuRE45["HepG2.2"].sum()



SuRE = pd.concat([SuRE42, SuRE43, SuRE44, SuRE45]).reset_index(drop=True)

SuRE["K562"] = np.mean(SuRE[["K562.1","K562.2","K562.3"]], axis=1)
SuRE["HepG2"] = np.mean(SuRE[["HepG2.1","HepG2.2"]], axis=1)

for individual in individuals:
    g = sns.lmplot(x="HepG2", y="K562", hue="rep", data=SuRE[SuRE["individual"] == individual], col="nodeClass")
    g.savefig(f"{individual}.pdf")


g = sns.lmplot(x="HepG2", y="K562", hue="rep", data=SuRE, col="nodeClass", row="individual")
g.savefig(f"regression.pdf")


SuRE = SuRE.groupby(["individual", "nodeClass", "name.ARBS"]).mean().reset_index()[["individual", "nodeClass", "name.ARBS", "lib", "K562.1","K562.2","K562.3","HepG2.1","HepG2.2", "K562", "HepG2"]]

"""
Note that we cannot do batch normalization, since these are not covers the whole batch,
  instead only ARBS region.
"""
# def batchNormalization(h0, e=0.001):
#   μ = h0.mean()
#   σ = h0.var()
#   ϵ = e
#   return (h0 - μ) / (σ - ϵ)**0.5
#
#
# SuRE["lib.1.n"] = batchNormalization(SuRE["lib"])
# SuRE["K562.1.n"] =  batchNormalization(SuRE["K562.1"])
# SuRE["K562.2.n"] =  batchNormalization(SuRE["K562.2"])
# SuRE["K562.3.n"] =  batchNormalization(SuRE["K562.3"])
# SuRE["HepG2.1.n"] =  batchNormalization(SuRE["HepG2.1"])
# SuRE["HepG2.2.n"] =  batchNormalization(SuRE["HepG2.2"])

SuRE[["K562.1.lib","K562.2.lib","K562.3.lib","HepG2.1.lib","HepG2.2.lib"]] = np.log10((SuRE[["K562.1","K562.2","K562.3","HepG2.1","HepG2.2"]].values + 1)/ SuRE["lib"].values[:,None])


SuRE[["K562.lib","HepG2.lib"]] = np.log10((SuRE[["K562","HepG2"]].values + 1)/ SuRE["lib"].values[:,None])

# SuRE42[["K562.1.lib","K562.2.lib","K562.3.lib","HepG2.1.lib","HepG2.2.lib"]] = (SuRE42[["K562.1","K562.2","K562.3","HepG2.1","HepG2.2"]].values + 1)/ SuRE42["lib"].values[:,None]
#
#
#
# SuRE43[["K562.1.lib","K562.2.lib","K562.3.lib","HepG2.1.lib","HepG2.2.lib"]] = (SuRE43[["K562.1","K562.2","K562.3","HepG2.1","HepG2.2"]].values + 1)/ SuRE43["lib"].values[:,None]
#
#
# SuRE44[["K562.1.lib","K562.2.lib","K562.3.lib","HepG2.1.lib","HepG2.2.lib"]] = (SuRE44[["K562.1","K562.2","K562.3","HepG2.1","HepG2.2"]].values + 1)/ SuRE44["lib"].values[:,None]
#
#
#
#
# SuRE45[["K562.1.lib","K562.2.lib","K562.3.lib","HepG2.1.lib","HepG2.2.lib"]] = (SuRE45[["K562.1","K562.2","K562.3","HepG2.1","HepG2.2"]].values + 1)/ SuRE45["lib"].values[:,None]
#
#




samples = ["SuRE42", "SuRE43", "SuRE44", "SuRE45"]
samples = [SuRE42, SuRE43, SuRE44, SuRE45]
#[(c,s) for s in samples for c in cols]



#[((x, hue), (x, hue))]



boxPairs = [("nAR", "con"),
            ("con", "ind"),
            ("con", "non"),
            ("ind", "non"),
            ("con", "pro")
        ]



colorPaletteAR = {"con": "#5A5A5A",
                  "ind": "#F9746D",
                  "non": "#ACACAC",
                  "nAR": "#D5D5D5",
                  "pro": "#000000"}



# cols = ["K562.1.lib","K562.2.lib","K562.3.lib","HepG2.1.lib","HepG2.2.lib"]
# individuals = ["SuRE42.1", "SuRE42.2", "SuRE43.1", "SuRE43.2", "SuRE44.1", "SuRE44.2", "SuRE45.1", "SuRE45.2"]
individuals = ["SuRE42", "SuRE43", "SuRE44", "SuRE45"]
cols = ["K562.lib","HepG2.lib"]
fig = plt.figure(figsize=[15, 10])
gs = gridspec.GridSpec(nrows=2, ncols=4)
plt.subplots_adjust(wspace=0.4, hspace=0.4)

params = dict(
    x="nodeClass",
    order=["con", "ind", "non", "nAR", "pro"]
)


i = 0
for individual in individuals:
    for c in cols:
        ax = fig.add_subplot(gs[i])
        ax = sns.boxplot(data=SuRE[SuRE["individual"] == individual],**params, y=c, palette=colorPaletteAR)
        m = int(np.ceil(SuRE[SuRE["individual"] == individual][c].max()))
        plt.title(f"{individual}\n{c}")
        plt.xlabel("NodeClasses", fontsize=14)
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(13)
        plt.ylabel("Sure Signal (log scale)", fontsize=14)
        ax.yaxis.set_ticks(range(m + 1))
        ax.yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
        ax.yaxis.set_ticks([np.log10(x) for p in range(m) for x in np.linspace(10**p, 10**(p+1), 10)], minor=True)
        add_stat_annotation(ax, data=SuRE[SuRE["individual"] == individual], **params, y=c, box_pairs=boxPairs, test='Mann-Whitney', loc='inside' ,line_height=0, text_offset=1,line_offset=0.05)
        i += 1

fig.tight_layout()
fig.savefig(f"SureFullMerged.pdf")

#######################################################################3





# cols = ["K562.1.lib","K562.2.lib","K562.3.lib","HepG2.1.lib","HepG2.2.lib"]
# individuals = ["SuRE42.1", "SuRE42.2", "SuRE43.1", "SuRE43.2", "SuRE44.1", "SuRE44.2", "SuRE45.1", "SuRE45.2"]
individuals = ["SuRE42", "SuRE43", "SuRE44", "SuRE45"]
replicates = [1, 2]
cols = ["K562.lib","HepG2.lib"]
fig = plt.figure(figsize=[15, 20])
gs = gridspec.GridSpec(nrows=4, ncols=4)
plt.subplots_adjust(wspace=0.4, hspace=0.4)

params = dict(
    x="nodeClass",
    order=["con", "ind", "non", "nAR", "pro"]
)


i = 0
for individual in individuals:
    for rep in replicates:
        f"{individual}.{rep}"
        for c in cols:
            ax = fig.add_subplot(gs[i])
            ax = sns.boxplot(data=SuRE[SuRE["individual"] == individual],**params, y=c, palette=colorPaletteAR)
            m = int(np.ceil(SuRE[SuRE["individual"] == individual][c].max()))
            plt.title(f"{individual}.{rep}\n{c}")
            plt.xlabel("NodeClasses", fontsize=14)
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(13)
            plt.ylabel("Sure Signal (log scale)", fontsize=14)
            ax.yaxis.set_ticks(range(m + 1))
            ax.yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
            ax.yaxis.set_ticks([np.log10(x) for p in range(m) for x in np.linspace(10**p, 10**(p+1), 10)], minor=True)
            add_stat_annotation(ax, data=SuRE[(SuRE["individual"] == individual) & (SuRE["rep"] == str(rep))], **params, y=c, box_pairs=boxPairs, test='Mann-Whitney', loc='inside' ,line_height=0, text_offset=1,line_offset=0.05)
            i += 1

fig.tight_layout()
fig.savefig(f"SureLibMerged.pdf")




############################3



home="/groups/lackgrp/ll_members/berkay/enhancerPromoterModel"

header = ["chr", "start", "end", "strand", "lib", "K562.B1.P1", "K562.B1.P2", "K562.B1.P3", "K562.B1.P4", "K562.B2.P1", "K562.B2.P2", "p.lib.1", "p.lib.2", "chr.ARBS", "start.ARBS", "end.ARBS", "name.ARBS"]

con = pd.read_table(home +"/GSE78709_SuRE-counts.con.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
con["nodeClass"] = "con"
ind = pd.read_table(home +"/GSE78709_SuRE-counts.ind.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
ind["nodeClass"] = "ind"
non = pd.read_table(home +"/GSE78709_SuRE-counts.non.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
non["nodeClass"] = "non"
nAR = pd.read_table(home +"/GSE78709_SuRE-counts.nAR.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
nAR["nodeClass"] = "nAR"


header = ["chr", "start", "end", "strand", "lib", "K562.B1.P1", "K562.B1.P2", "K562.B1.P3", "K562.B1.P4", "K562.B2.P1", "K562.B2.P2", "p.lib.1", "p.lib.2", "chr.ARBS", "start.ARBS", "end.ARBS", "name.ARBS", "strand.ARBS"]

pro = pd.read_table(home +"/GSE78709_SuRE-counts.tss.bed", names=header).groupby(["name.ARBS"]).sum().reset_index()
pro["nodeClass"] = "pro"

SuRE1 = pd.concat([con, ind, non, nAR, pro])

SuRE1[["K562.B1.P1.lib", "K562.B1.P2.lib", "K562.B1.P3.lib", "K562.B1.P4.lib", "K562.B2.P1.lib", "K562.B2.P2.lib"]] = np.log10((SuRE1[["K562.B1.P1", "K562.B1.P2", "K562.B1.P3", "K562.B1.P4", "K562.B2.P1", "K562.B2.P2"]].values + 1)/ SuRE1["lib"].values[:,None])




boxPairs = [("nAR", "con"),
            ("con", "ind"),
            ("con", "non"),
            ("ind", "non"),
            ("con", "pro")
        ]



colorPaletteAR = {"con": "#5A5A5A",
                  "ind": "#F9746D",
                  "non": "#ACACAC",
                  "nAR": "#D5D5D5",
                  "pro": "#000000"}



cols = ["K562.B1.P1.lib", "K562.B1.P2.lib", "K562.B1.P3.lib", "K562.B1.P4.lib", "K562.B2.P1.lib", "K562.B2.P2.lib"]

fig = plt.figure(figsize=[8,8])
gs = gridspec.GridSpec(ncols=3, nrows=2)
plt.subplots_adjust(wspace=0.4, hspace=0.4)

params = dict(
    x="nodeClass",
    order=["con", "ind", "non", "nAR", "pro"]
)


i = 0
for c in cols:
    ax = fig.add_subplot(gs[i])
    ax = sns.boxplot(data=SuRE1,**params, y=c, palette=colorPaletteAR)
    m = int(np.ceil(SuRE1[c].max()) + 1)
    plt.title(c)
    plt.xlabel("NodeClasses", fontsize=14)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(13)
    plt.ylabel("Sure Signal (log scale)", fontsize=14)
    ax.yaxis.set_ticks(range(m))
    ax.yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
    ax.yaxis.set_ticks([np.log10(x) for p in range(m) for x in np.linspace(10**p, 10**(p+1), 10)], minor=True)
    add_stat_annotation(ax, data=SuRE1, **params, y=c, box_pairs=boxPairs, test='Mann-Whitney', loc='inside' ,line_height=0, text_offset=1,line_offset=0.05)
    i += 1


fig.savefig(f"Sure1.png")
