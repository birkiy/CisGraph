

from Functions.Packages import *

home = "/home/ualtintas/GR/histoneSignal"
header = ["H3K4me1.Etoh","H3K4me1.Dex","H3K27ac.Etoh","H3K27ac.Dex", "H3K4me3.Etoh", "H3K4me3.Dex", "RNAP2.Etoh", "RNAP2.Dex"]
con = pd.read_csv(f"{home}/GR.con.tab", sep="\t", names=header)
con["nodeClass"] = "con"
ind = pd.read_csv(f"{home}/GR.ind.tab", sep="\t", names=header)
ind["nodeClass"] = "ind"
nAR = pd.read_csv(f"{home}/negativeControlGR.final.tab", sep="\t", names=header)
nAR["nodeClass"] = "nAR"
non = pd.read_csv(f"{home}/GR.non.tab", sep="\t", names=header)
non["nodeClass"] = "non"
tss = pd.read_csv(f"{home}/TSS.hg19.Idx.woS.tab", sep="\t", names=header)
tss["nodeClass"] = "tss"


DF = pd.concat([con, ind, non, nAR, tss])


NDF = sklearn.preprocessing.normalize(DF[header])
NDF = pd.DataFrame(NDF, columns=header, index=DF.index)
NDF["nodeClass"] = DF["nodeClass"]

fig = plt.figure(figsize=[8, 5])
gs = gridspec.GridSpec(ncols=2, nrows=1)
plt.subplots_adjust(wspace=0.4)


colorPalette = ["#63b7af", "#abf0e9","#d4f3ef", "#f5fffd", "#ee8572"]

params = dict(x='nodeClass',
              order=["con", "ind", "non", "nAR", "tss"],
              data=NDF)


boxPairs = [("nAR", "con"),
            ("con", "ind"),
            ("con", "non"),
            ("ind", "non"),
            ("con", "tss"),
            ("non", "tss"),
            ("nAR", "tss")
            ]

target = "H3K27ac"
ax1 = fig.add_subplot(gs[0])
ax1 = sns.boxplot(**params, y=f"{target}.Etoh", palette=colorPalette)
plt.title(f"{target}.Etoh")
plt.ylabel(f"{target} Signal Level", fontsize=14)
plt.xlabel("Node Classes", fontsize=14)
add_stat_annotation(ax1, **params, y=f"{target}.Etoh", box_pairs=boxPairs, test='Mann-Whitney')

ax2 = fig.add_subplot(gs[1])
ax2 = sns.boxplot(**params, y=f"{target}.Dex", palette=colorPalette)
plt.title(f"{target}.Dex")
plt.ylabel(f"{target} Signal Level", fontsize=14)
plt.xlabel("Node Classes", fontsize=14)
add_stat_annotation(ax2, **params, y=f"{target}.Dex", box_pairs=boxPairs, test='Mann-Whitney')

fig.savefig(f"{figureRoot}/GR.{target}.pdf")
