
files = [
"A549.H3K27ac.dex.0h.rep1.con.avg.tab",
"A549.H3K27ac.dex.0h.rep1.ind.avg.tab",
"A549.H3K27ac.dex.0h.rep1.non.avg.tab",
"A549.H3K27ac.dex.0h.rep2.con.avg.tab",
"A549.H3K27ac.dex.0h.rep2.ind.avg.tab",
"A549.H3K27ac.dex.0h.rep2.non.avg.tab",
"A549.H3K27ac.dex.0h.rep3.ind.avg.tab",
"A549.H3K27ac.dex.8h.rep1.con.avg.tab",
"A549.H3K27ac.dex.8h.rep1.ind.avg.tab",
"A549.H3K27ac.dex.8h.rep1.non.avg.tab",
"A549.H3K27ac.dex.8h.rep2.con.avg.tab",
"A549.H3K27ac.dex.8h.rep2.ind.avg.tab",
"A549.H3K27ac.dex.8h.rep2.non.avg.tab",
"A549.H3K4me1.dex.0h.rep1.ind.avg.tab",
"A549.H3K4me1.dex.0h.rep2.ind.avg.tab",
"A549.H3K4me1.dex.0h.rep3.con.avg.tab",
"A549.H3K4me1.dex.0h.rep3.ind.avg.tab",
"A549.H3K4me1.dex.0h.rep3.non.avg.tab",
"A549.H3K4me1.dex.8h.rep1.ind.avg.tab",
"A549.H3K4me1.dex.8h.rep2.con.avg.tab",
"A549.H3K4me1.dex.8h.rep2.ind.avg.tab",
"A549.H3K4me1.dex.8h.rep2.non.avg.tab",
"A549.H3K4me1.dex.8h.rep3.ind.avg.tab",
"A549.H3K4me3.dex.0h.rep1.con.avg.tab",
"A549.H3K4me3.dex.0h.rep1.ind.avg.tab",
"A549.H3K4me3.dex.0h.rep1.non.avg.tab",
"A549.H3K4me3.dex.0h.rep2.con.avg.tab",
"A549.H3K4me3.dex.0h.rep2.ind.avg.tab",
"A549.H3K4me3.dex.0h.rep2.non.avg.tab",
"A549.H3K4me3.dex.0h.rep3.ind.avg.tab",
"A549.H3K4me3.dex.8h.rep1.ind.avg.tab",
"A549.H3K4me3.dex.8h.rep2.con.avg.tab",
"A549.H3K4me3.dex.8h.rep2.ind.avg.tab",
"A549.H3K4me3.dex.8h.rep2.non.avg.tab",
"A549.POLR2A.dex.1h.rep1.con.avg.tab",
"A549.POLR2A.dex.1h.rep1.ind.avg.tab",
"A549.POLR2A.dex.1h.rep1.non.avg.tab",
"A549.POLR2A.dex.1h.rep2.con.avg.tab",
"A549.POLR2A.dex.1h.rep2.ind.avg.tab",
"A549.POLR2A.dex.1h.rep2.non.avg.tab",
"A549.POLR2A.etoh.1h.rep1.con.avg.tab",
"A549.POLR2A.etoh.1h.rep1.ind.avg.tab",
"A549.POLR2A.etoh.1h.rep1.non.avg.tab",
"A549.POLR2A.etoh.1h.rep2.con.avg.tab",
"A549.POLR2A.etoh.1h.rep2.ind.avg.tab",
"A549.POLR2A.etoh.1h.rep2.non.avg.tab"
]




import pandas as pd
import seaborn as sns
from statannot import add_stat_annotation
import matplotlib.pyplot as plt
from matplotlib import gridspec, colors
import numpy as np
from matplotlib import ticker as mticker


home = "/groups/lackgrp/ll_members/berkay/STARRbegin/A549"
file = f"{home}/coverage/A549.H3K27ac.dex.0h.rep1.con.avg.tab"


def unPosName(DF):
    new = DF["name"].str.split(":", expand=True)
    DF["chr"] = new[0].str.strip()
    new = new[1].str.split("-", expand=True)
    DF["start"] = new[0].str.strip().astype("int64")
    DF["end"] = new[1].str.strip().astype("int64")
    return DF

def processBAOB(file, raw):
    header = ["name", "size", "covered", "sum", "mean0", "mean"]
    DF = pd.read_table(file, names=header)
    DF = unPosName(DF)
    new = raw.split(".")
    DF["target"] = new[1]
    DF["treatment"] = new[2]
    DF["duration"] = new[3]
    DF["rep"] = new[4]
    DF["nodeClass"] = new[5]
    return DF[["chr", "start", "end", "name", "nodeClass", "target", "treatment", "duration", "rep", "sum", "mean"]]



DF = pd.DataFrame()
for file in files:
    DF = pd.concat([DF, processBAOB(f"{home}/coverage/{file}", file)])


DF = DF.reset_index(drop=True)

DFpivot = np.log(DF.pivot_table(values="mean", columns="target", index=["nodeClass", "treatment", "duration", "name"]) +1 ).reset_index().fillna(0)





rows = ["0h", "8h"]
cols = ["H3K4me3","H3K27ac","H3K4me1"]
fig = plt.figure(figsize=[8, 8])
gs = gridspec.GridSpec(nrows=3, ncols=2)
plt.subplots_adjust(wspace=0.4, hspace=0.4)

params = dict(
    x="nodeClass",
    order=["con", "ind", "non", "nAR", "pro"]
)

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



i = 0
for r in rows:
    for c in cols:
        ax = fig.add_subplot(gs[i])
        ax = sns.boxplot(data=DFpivot[DFpivot["duration"] == r], y=c,**params, palette=colorPaletteAR)
        m = int(np.ceil(DFpivot[DFpivot["duration"] == r][c].max()))
        plt.title(f"{r}\n{c}")
        plt.xlabel("NodeClasses", fontsize=14)
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(13)
        plt.ylabel("Occupancy Score (log scale)", fontsize=14)
        ax.yaxis.set_ticks(range(m + 1))
        ax.yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
        ax.yaxis.set_ticks([np.log10(x) for p in range(m) for x in np.linspace(10**p, 10**(p+1), 10)], minor=True)
        add_stat_annotation(ax, data=DFpivot[DFpivot["duration"] == r], **params, y=c, box_pairs=boxPairs, test='Mann-Whitney', loc='inside' ,line_height=0, text_offset=1,line_offset=0.05)
        print(r, c)
        i += 1

fig.tight_layout()
fig.savefig(f"Histones.pdf")
