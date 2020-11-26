
from Functions.Packages import *


GroDF = pickle.load(open(f"{dataRoot}/tmpData/Gro.DF.p", "rb" ))
GroDF["logFC"] = None
GroDF["Qval"] = None
GroDF = GroDF.reset_index()
GroDF = GroDF.rename(columns={"index": "Gene"})
GroDF = GroDF.sort_values("Gene")
GroDF = GroDF.reset_index().drop("index", axis=1)


logFC = pd.read_csv(f"{dataRoot}/DEG/GSE64529_diffexpr-results.csv")
logFC = logFC[["Gene", "log2FoldChange", "padj"]]
logFC = logFC.sort_values("Gene").reset_index().drop("index", axis=1)
logFC["Gene"] = logFC["Gene"].astype(str)

i = 0
for gene, lfc, q in zip(logFC["Gene"], logFC["log2FoldChange"], logFC["padj"]):
    left = bi.bisect_left(GroDF["Gene"] , gene +".")
    right = bi.bisect_left(GroDF["Gene"] , gene+".z")
    GroDF.loc[left:right-1, "logFC"] = lfc
    GroDF.loc[left:right-1, "Qval"] = q
    i += 1
    if i % 1000 == 0:
        print(gene)
        print(GroDF.loc[left:right-1, "Gene"])
        print(i)




conAnn = pd.read_csv("/home/ualtintas/ARBSs/annotations/con.Annotations.txt", sep="\t")
conAnn = conAnn.rename(columns={"PeakID (cmd=annotatePeaks.pl regions/cons-arbs.bed hg19)": "PeakID"})
indAnn = pd.read_csv("/home/ualtintas/ARBSs/annotations/ind.Annotations.txt", sep="\t")
indAnn = indAnn.rename(columns={"PeakID (cmd=annotatePeaks.pl regions/ind-arbs.bed hg19)": "PeakID"})
nonAnn = pd.read_csv("/home/ualtintas/ARBSs/annotations/non.Annotations.txt", sep="\t")
nonAnn = nonAnn.rename(columns={"PeakID (cmd=annotatePeaks.pl regions/Non-Active-ARBS.bed hg19)": "PeakID"})





capG = pd.concat([con, ind, non]).reset_index(drop=True)

capG["annotation"] = capG["annotation"].str.split("(", expand=True)[0]


GroARBS = GroDF.loc[GroDF["Gene"].isin(capG["PeakID"]) , ["Gene","dmso.-", "dmso.+", "dht.-", "dht.+"]].reset_index(drop=True)
GroARBS = GroARBS.rename(columns={"Gene": "PeakID"})
capG = capG.sort_values("PeakID").reset_index(drop=True)

if sum(GroARBS["PeakID"] == capG["PeakID"]) == capG.shape[0]:
    signals = GroARBS.iloc[:,1:]
    capG = pd.concat([capG,signals], axis=1)



capG["NumberofCapsDHT.-"] += 1
capG["NumberofCapsDHT.+"] += 1

capG["NumberofCapsVEH.-"] += 1
capG["NumberofCapsVEH.+"] += 1



fig = plt.figure(figsize=[15,15])
plt.subplots_adjust(bottom = 0.2)
gs = gridspec.GridSpec(ncols=2, nrows=2)



fig.add_subplot(gs[0,0])
sns.scatterplot(x="SignalVEH.+", y="dmso.+", data=capG, hue="nodeClass", size="NumberofCapsVEH.+", hue_order=["con", "ind", "non"])
plt.legend()
# plt.yscale("log")
# plt.xscale("log")

fig.add_subplot(gs[0,1])
sns.scatterplot(x="SignalDHT.+", y="dht.+", data=capG, hue="nodeClass", size="NumberofCapsDHT.+", hue_order=["con", "ind", "non"])
plt.legend()
# plt.yscale("log")
# plt.xscale("log")

fig.add_subplot(gs[1,0])
sns.scatterplot(x="SignalVEH.-", y="dmso.-", data=capG, hue="nodeClass", size="NumberofCapsVEH.-", hue_order=["con", "ind", "non"])
plt.legend()
# plt.yscale("log")
# plt.xscale("log")

fig.add_subplot(gs[1,1])
sns.scatterplot(x="SignalDHT.-", y="dht.-", data=capG, hue="nodeClass", size="NumberofCapsDHT.-", hue_order=["con", "ind", "non"])
plt.legend()
# plt.yscale("log")
# plt.xscale("log")

fig.savefig(f"{figureRoot}/ProCapvsGro.ARBS.pdf")



df2 = capG.groupby(["nodeClass", "Annotation"]).size().reset_index()

df2.pivot(index="Annotation", columns="nodeClass", values=0).fillna(0)
