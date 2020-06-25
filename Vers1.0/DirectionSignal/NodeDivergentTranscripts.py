
from Functions.Packages import *
from Functions.PlotFunctions import *

G = pickle.load(open(f"{dataRoot}/tmpData/GraphsGData.p", "rb" ))


nodeDTlvl = pd.DataFrame()
for i, node in enumerate(G.nodes()):
    if G.nodes[node]["nodeClass"] in ("pro", "dwP", "upP"):
        continue
    PlvlDms = G.nodes[node]["lvlP.DMSO"]
    MlvlDms = G.nodes[node]["lvlM.DMSO"]
    PlvlDht = G.nodes[node]["lvlP.DHT"]
    MlvlDht = G.nodes[node]["lvlM.DHT"]
    nodeClass = G.nodes[node]["nodeClass"]
    nodeDTlvl = pd.concat([nodeDTlvl, pd.DataFrame({"lvlP.DMSO": PlvlDms,"lvlP.DHT": PlvlDht, "lvlM.DMSO": MlvlDms,"lvlM.DHT": MlvlDht, "nodeClass": nodeClass}, index=[i])])





nodeDTlvl["logPlvl.DMSO"] = np.log(nodeDTlvl["lvlP.DMSO"] )
nodeDTlvl["logPlvl.DHT"] = np.log(nodeDTlvl["lvlP.DHT"] )
nodeDTlvl["logMlvl.DMSO"] = np.log(nodeDTlvl["lvlM.DMSO"] )
nodeDTlvl["logMlvl.DHT"] = np.log(nodeDTlvl["lvlM.DHT"] )

nodeDTlvl.loc[nodeDTlvl["nodeClass"].isin(["tsP", "uPP", "dPP"]), "nodeClass"] = "tsP"
nodeDTlvl.loc[nodeDTlvl["nodeClass"].isin(["tsM", "uMP", "dMP"]), "nodeClass"] = "tsM"

nodeDTlvl = nodeDTlvl[~(nodeDTlvl["nodeClass"].isin([ "pro", "upP", "dwP", "tsP", "tsM"]))]


filtDF = nodeDTlvl[~((nodeDTlvl["logPlvl.DMSO"] == 0) & (nodeDTlvl["logMlvl.DMSO"] == 0) & (nodeDTlvl["logPlvl.DHT"] == 0) & (nodeDTlvl["logMlvl.DHT"] == 0))]

dmsDF = filtDF[["lvlP.DMSO", "lvlM.DMSO", "logPlvl.DMSO", "logMlvl.DMSO", "nodeClass"]]
dmsDF = dmsDF.rename(columns={"lvlP.DMSO": "lvlP", "lvlM.DMSO": "lvlM", "logPlvl.DMSO": "logPlvl", "logMlvl.DMSO": "logMlvl"})
dmsDF["Condition"] = "DMSO"
dhtDF = filtDF[["lvlP.DHT", "lvlM.DHT", "logPlvl.DHT", "logMlvl.DHT", "nodeClass"]]
dhtDF = dhtDF.rename(columns={"lvlP.DHT": "lvlP", "lvlM.DHT": "lvlM", "logPlvl.DHT": "logPlvl", "logMlvl.DHT": "logMlvl"})
dhtDF["Condition"] = "DHT"


DF = pd.concat([dmsDF, dhtDF])




order=["con", "ind", "non", "oth", "uPP", "dPP", "uMP", "dMP"]
hueOrder = ["DMSO", "DHT"]
colorPalette = {"DMSO": "#FAB550", "DHT": "#5F9BD2"}

combs = [
    (a, b)
    for a in order
    for b in hueOrder
    ]

boxPairs = [
    (a, b)
    for a in combs
    for b in combs
    if a != b
]

fig = plt.figure(figsize=[10,5])
gs = gridspec.GridSpec(ncols=2, nrows=1)

fig.add_subplot(gs[0])
ax1 = sns.boxplot(data=DF, x="nodeClass", y="logPlvl", hue="Condition", order=order, hue_order=hueOrder, palette=colorPalette)
# add_stat_annotation(ax1, data=DF, x="nodeClass", y="logPlvl", hue="Condition", box_pairs=boxPairs, test='Mann-Whitney', loc='inside')

fig.add_subplot(gs[1])
ax2 = sns.boxplot(data=DF, x="nodeClass", y="logMlvl", hue="Condition", order=order, hue_order=hueOrder, palette=colorPalette)
# add_stat_annotation(ax2, data=DF, x="nodeClass", y="logMlvl", hue="Condition", box_pairs=boxPairs, test='Mann-Whitney', loc='inside')

fig.savefig(f"{figureRoot}/nodeTlvls.filter.pdf")


plt.close("all")

fig = plt.figure(figsize=[12,9])
gs = gridspec.GridSpec(ncols=2, nrows=1)

fig.add_subplot(gs[0])
ax1 = sns.boxplot(data=DF, x="Condition", y="logPlvl", hue="nodeClass", order=hueOrder, hue_order=order)
# add_stat_annotation(ax1, data=DF, x="nodeClass", y="logPlvl", hue="Condition", box_pairs=boxPairs, test='Mann-Whitney', loc='inside')

fig.add_subplot(gs[1])
ax2 = sns.boxplot(data=DF, x="Condition", y="logMlvl", hue="nodeClass", order=hueOrder, hue_order=order)
# add_stat_annotation(ax2, data=DF, x="nodeClass", y="logMlvl", hue="Condition", box_pairs=boxPairs, test='Mann-Whitney', loc='inside')

fig.savefig(f"{figureRoot}/nodeTlvls.filter.nodeHue.pdf")



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



filtGroDF = GroDF[~(
((GroDF["logFC"] < 1) | (GroDF["logFC"] > -1)) &
(GroDF["Qval"] > 0.05) &
(GroDF["class"].isin(["tsP", "tsM"]))
)].sort_values("logFC").reset_index().drop("index", axis=1)


filtGroDF.loc[((filtGroDF["logFC"] > 1) & (filtGroDF["class"] == "tsP")), "class"] = "uPP"


filtGroDF.loc[((filtGroDF["logFC"] > 1) & (filtGroDF["class"] == "tsM")), "class"] = "uMP"


filtGroDF.loc[((filtGroDF["logFC"] < -1) & (filtGroDF["class"] == "tsP")), "class"] = "dPP"


filtGroDF.loc[((filtGroDF["logFC"] < -1) & (filtGroDF["class"] == "tsM")), "class"] = "dMP"




dmsDF = filtGroDF[["dmso.-", "dmso.+", "class"]]
dmsDF = dmsDF.rename(columns={"dmso.-": "GroSeqMinus", "dmso.+": "GroSeqPlus"})
dmsDF["Condition"] = "DMSO"
dhtDF = filtGroDF[["dht.-", "dht.+", "class"]]
dhtDF = dhtDF.rename(columns={"dht.-": "GroSeqMinus", "dht.+": "GroSeqPlus"})
dhtDF["Condition"] = "DHT"


allDF = pd.concat([dmsDF, dhtDF])


hueOrder = ["DMSO", "DHT"]
colorPalette = {"DMSO": "#FAB550", "DHT": "#5F9BD2"}
order=["con", "ind", "non", "oth", "uPP", "dPP", "uMP", "dMP"]


fig = plt.figure(figsize=[10,5])
gs = gridspec.GridSpec(ncols=2, nrows=1)

fig.add_subplot(gs[0])
ax1 = sns.boxplot(data=allDF, x="class", y="GroSeqPlus", hue="Condition", order=order, hue_order=hueOrder, palette=colorPalette)
# add_stat_annotation(ax1, data=DF, x="nodeClass", y="logPlvl", hue="Condition", box_pairs=boxPairs, test='Mann-Whitney', loc='inside')
ax1.set(yscale="log")

fig.add_subplot(gs[1])
ax2 = sns.boxplot(data=allDF, x="class", y="GroSeqMinus", hue="Condition", order=order, hue_order=hueOrder, palette=colorPalette)
# add_stat_annotation(ax2, data=DF, x="nodeClass", y="logMlvl", hue="Condition", box_pairs=boxPairs, test='Mann-Whitney', loc='inside')
ax2.set(yscale="log")

fig.savefig(f"{figureRoot}/nodeTlvls.allDE.pdf")


plt.close("all")



























geneNodes = list(GroDF[GroDF["class"].isin(["tsP", "tsM"])].index)


fileCSV = f"{dataRoot}/DEG/GSE64529_diffexpr-results.csv"
with open(fileCSV) as csvFile:
    reader = csv.reader(csvFile, delimiter=',')
    next(reader)
    for i, row in enumerate(reader):
        geneSymbol = row[1]
        logFC = row[3]
        Qval = row[7]
        if (logFC == "NA" or Qval == "NA"):
            continue
        logFC = float(logFC)
        Qval = float(Qval)
        GroDF.loc[GroDF["index"].str.find(geneSymbol + ".") != -1, "logFC"] = logFC
        GroDF.loc[GroDF["index"].str.find(geneSymbol + ".") != -1, "Qval"] = Qval
        if i % 1000 == 0:
            print(i)




filtGroDF = GroDF[~((GroDF["logPlvl.DMSO"] == 0) & (GroDF["logMlvl.DMSO"] == 0) & (GroDF["logPlvl.DHT"] == 0) & (GroDF["logMlvl.DHT"] == 0))]

dmsDF = GroDF[["dmso.-", "dmso.+", "class"]]
dmsDF = dmsDF.rename(columns={"dmso.-": "GroSeqMinus", "dmso.+": "GroSeqPlus"})
dmsDF["Condition"] = "DMSO"
dhtDF = GroDF[["dht.-", "dht.+", "class"]]
dhtDF = dhtDF.rename(columns={"dht.-": "GroSeqMinus", "dht.+": "GroSeqPlus"})
dhtDF["Condition"] = "DHT"


allDF = pd.concat([dmsDF, dhtDF])


order=["con", "ind", "non", "oth", "tsP", "tsM"]


fig = plt.figure(figsize=[10,5])
gs = gridspec.GridSpec(ncols=2, nrows=1)

fig.add_subplot(gs[0])
ax1 = sns.boxplot(data=allDF, x="class", y="GroSeqPlus", hue="Condition", order=order, hue_order=hueOrder, palette=colorPalette)
# add_stat_annotation(ax1, data=DF, x="nodeClass", y="logPlvl", hue="Condition", box_pairs=boxPairs, test='Mann-Whitney', loc='inside')
ax1.set(yscale="log")

fig.add_subplot(gs[1])
ax2 = sns.boxplot(data=allDF, x="class", y="GroSeqMinus", hue="Condition", order=order, hue_order=hueOrder, palette=colorPalette)
# add_stat_annotation(ax2, data=DF, x="nodeClass", y="logMlvl", hue="Condition", box_pairs=boxPairs, test='Mann-Whitney', loc='inside')
ax2.set(yscale="log")

fig.savefig(f"{figureRoot}/nodeTlvls.all.pdf")


plt.close("all")
