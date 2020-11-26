

from Functions.Packages import *


GroDF = pickle.load(open(f"{dataRoot}/Gro/Gro.DF.p", "rb" ))
GroDF["logFC"] = None
GroDF["Qval"] = None
GroDF = GroDF.reset_index()
GroDF = GroDF.rename(columns={"index": "Gene"})
GroDF = GroDF.sort_values("Gene")
GroDF = GroDF.reset_index().drop("index", axis=1)

proBed = readBed(f"/home/ualtintas/genomeAnnotations/Regions/TSS.hg19.Idx.bed")
proBed = pd.DataFrame.from_dict(proBed, orient="index").reset_index()
proBed = proBed.rename(columns={"index": "Gene", 0: "chr", 1: "start", 2: "end"})
proBed = proBed.sort_values("Gene").reset_index(drop=True)
proBed["logFC"] = None
proBed["Qval"] = None


logFC = pd.read_csv(f"{dataRoot}/DEG/GSE64529_diffexpr-results.csv")
logFC = logFC[["Gene", "log2FoldChange", "padj"]]
logFC = logFC.sort_values("Gene").reset_index().drop("index", axis=1)
logFC["Gene"] = logFC["Gene"].astype(str)

i = 0
for gene, lfc, q in zip(logFC["Gene"], logFC["log2FoldChange"], logFC["padj"]):
    left = bi.bisect_left(proBed["Gene"] , gene +".")
    right = bi.bisect_left(proBed["Gene"] , gene+".z")
    proBed.loc[left:right-1, "logFC"] = lfc
    proBed.loc[left:right-1, "Qval"] = q
    i += 1
    if i % 1000 == 0:
        print(gene)
        print(proBed.loc[left:right-1, "Gene"])
        print(i)


proDE = proBed[(((proBed["logFC"] > 1) | (proBed["logFC"] < -1)) & (proBed["Qval"] < 0.05))]

proDE = proDE[["chr", "start", "end", "Gene"]]

proDE.to_csv(f"{dataRoot}/Regions/proDE.bed", sep="\t", index=False, header=False)

conAnn = pd.read_csv(f"{dataRoot}/Regions/con.AnnFile.cvs", sep="\t")
indAnn = pd.read_csv(f"{dataRoot}/Regions/ind.AnnFile.csv", sep="\t")
nonAnn = pd.read_csv(f"{dataRoot}/Regions/non.AnnFile.csv", sep="\t")
nonAnn = pd.read_csv(f"{dataRoot}/Regions/nAR.AnnFile.csv", sep="\t")
othAnn = pd.read_csv(f"{dataRoot}/Regions/oth.AnnFile.csv", sep="\t")
tsPAnn = pd.read_csv("/home/ualtintas/genomeAnnotations/Regions/TSS.hg19.+.AnnFile.csv", sep="\t")
tsMAnn = pd.read_csv("/home/ualtintas/genomeAnnotations/Regions/TSS.hg19.-.AnnFile.csv", sep="\t")

conAnn = conAnn.rename(columns={"seqnames": "Chr", "start": "Start", "end": "End"})
indAnn = indAnn.rename(columns={"seqnames": "Chr", "start": "Start", "end": "End"})
nonAnn = nonAnn.rename(columns={"seqnames": "Chr", "start": "Start", "end": "End"})
othAnn = othAnn.rename(columns={"seqnames": "Chr", "start": "Start", "end": "End"})
tsPAnn = tsPAnn.rename(columns={"seqnames": "Chr", "start": "Start", "end": "End"})
tsMAnn = tsMAnn.rename(columns={"seqnames": "Chr", "start": "Start", "end": "End"})



con = conAnn[["V4", "Chr", "Start", "End", "annotation", "geneStrand", "SYMBOL"]]
ind = indAnn[["V4", "Chr", "Start", "End", "annotation", "geneStrand", "SYMBOL"]]
non = nonAnn[["V4", "Chr", "Start", "End", "annotation", "geneStrand", "SYMBOL"]]
oth = othAnn[["V4", "Chr", "Start", "End", "annotation", "geneStrand", "SYMBOL"]]
tsP = tsPAnn[["V4", "Chr", "Start", "End", "annotation", "geneStrand", "SYMBOL"]]
tsM = tsMAnn[["V4", "Chr", "Start", "End", "annotation", "geneStrand", "SYMBOL"]]


con["nodeClass"] = "con"
ind["nodeClass"] = "ind"
non["nodeClass"] = "non"
oth["nodeClass"] = "oth"
tsP["nodeClass"] = "tsP"
tsM["nodeClass"] = "tsM"


# allAnn = pd.concat([con, ind, non, oth, tsP, tsM]). reset_index(drop=True)
allAnn = pd.concat([con, ind, non]). reset_index(drop=True)


allAnn["annotation"] = allAnn["annotation"].str.split("(", expand=True)[0]



i = 0
for peak, ann, st, sym in zip(allAnn["V4"], allAnn["annotation"], allAnn["geneStrand"], allAnn["SYMBOL"]):
    left = bi.bisect_left(GroDF["Gene"] , peak)
    right = bi.bisect_left(GroDF["Gene"] , peak+"z")
    GroDF.loc[left:right-1, "annotation"] = ann
    GroDF.loc[left:right-1, "geneStrand"] = st
    GroDF.loc[left:right-1, "SYMBOL"] = sym
    i += 1
    if i % 1000 == 0:
        print(peak)
        print(GroDF.loc[left:right-1, "Gene"])
        print(i)







dmsDF = GroDF[["Gene", "dmso.-", "dmso.+", "class", "annotation", "geneStrand"]]
dmsDF = dmsDF.rename(columns={"dmso.-": "GroSeqMinus", "dmso.+": "GroSeqPlus"})
dmsDF["Condition"] = "DMSO"
dhtDF = GroDF[["Gene", "dht.-", "dht.+", "class", "annotation", "geneStrand"]]
dhtDF = dhtDF.rename(columns={"dht.-": "GroSeqMinus", "dht.+": "GroSeqPlus"})
dhtDF["Condition"] = "DHT"




dmsDF = GroDF[["Gene", "dmso.-", "dmso.+", "class"]]
dmsDF = dmsDF.rename(columns={"dmso.-": "GroSeqMinus", "dmso.+": "GroSeqPlus"})
dmsDF["Condition"] = "DMSO"
dhtDF = GroDF[["Gene", "dht.-", "dht.+", "class"]]
dhtDF = dhtDF.rename(columns={"dht.-": "GroSeqMinus", "dht.+": "GroSeqPlus"})
dhtDF["Condition"] = "DHT"








df = dmsDF.groupby(["annotation", "class"]).size()

df.pivot(index="annotation", columns="class", values=0)



allDF = pd.concat([dmsDF, dhtDF])

# Final Final


allDF["GroSeqSum"] = allDF["GroSeqMinus"] + allDF["GroSeqMinus"]


allDF = allDF.rename(columns={"class": "nodeClass"})

allDF.loc[allDF["nodeClass"] == "tsP", "nodeClass"] = "pro"
allDF.loc[allDF["nodeClass"] == "tsM", "nodeClass"] = "pro"



fig = plt.figure(figsize=[9, 7])
gs = gridspec.GridSpec(ncols=2, nrows=1)
# plt.subplots_adjust(top = 0.8, hspace=0.5)


boxPairs = [("tss", "con"),
            ("tss", "ind"),
            ("tss", "non"),
            ("nAR", "non"),
            ("nAR", "ind"),
            ("nAR", "tss"),
            ("nAR", "con"),
            ("con", "ind"),
            ("con", "non"),
            ("ind", "non")
            ]


colorPalette = ["#63b7af", "#abf0e9","#d4f3ef", "#f5fffd", "#ee8572"]

params = dict(
              x='nodeClass',
              y='GroSeqSum',
              order=["con", "ind", "non", "nAR", "tss"])


ax1 = fig.add_subplot(gs[0])
ax1 = sns.boxplot(**params, data=allDF[allDF["Condition"] == "DMSO"], palette=colorPalette)
plt.title("DMSO")
ax1.set(yscale="log")


add_stat_annotation(ax1, **params, box_pairs=boxPairs, data=allDF[allDF["Condition"] == "DMSO"], test='Mann-Whitney')

ax2 = fig.add_subplot(gs[1])
ax2 = sns.boxplot(**params, data=allDF[allDF["Condition"] == "DHT"], palette=colorPalette)
plt.title("DHT")
ax2.set(yscale="log")

add_stat_annotation(ax2, **params, box_pairs=boxPairs, data=allDF[allDF["Condition"] == "DHT"], test='Mann-Whitney')


fig.savefig(f"{figureRoot}/Gro.ARBS.pdf")






















import itertools

order=["con", "ind", "non", "oth", "tsP", "tsM"]
hueOrder = ["DMSO", "DHT"]
boxPairs = list(itertools.combinations([(x, h) for x in order for h in hueOrder], 2))
colorPalette = {"DMSO": "#FAB550", "DHT": "#5F9BD2"}

fig = plt.figure(figsize=[10, 5])
gs = gridspec.GridSpec(ncols=2, nrows=1)

fig.add_subplot(gs[0])
ax1 = sns.boxplot(data=allDF, x="class", y="GroSeqPlus", hue="Condition", order=order, hue_order=hueOrder, palette=colorPalette)
# add_stat_annotation(ax1, data=allDF, x="class", y="GroSeqPlus", hue="Condition", box_pairs=boxPairs, test='Mann-Whitney', loc='outside')
ax1.set(yscale="log")

fig.add_subplot(gs[1])
ax2 = sns.boxplot(data=allDF, x="class", y="GroSeqMinus", hue="Condition", order=order, hue_order=hueOrder, palette=colorPalette)
# add_stat_annotation(ax2, data=allDF, x="class", y="GroSeqMinus", hue="Condition", box_pairs=boxPairs, test='Mann-Whitney', loc='outside')
ax2.set(yscale="log")

fig.savefig(f"{figureRoot}/Gro.all.pdf")

plt.close("all")





hueOrder=["con", "ind", "non", "oth", "tsP", "tsM"]
order = ["3' UTR", "5' UTR", "Distal Intergenic", "Downstream ",  "Exon ", "Intron ", "Promoter "]

fig = plt.figure(figsize=[20, 40])
gs = gridspec.GridSpec(ncols=2, nrows=4)
plt.subplots_adjust(bottom = 0.2, hspace=0.5)


fig.add_subplot(gs[0])
ax1 = sns.boxplot(data=dmsDF[dmsDF["geneStrand"] == 1], x="annotation", y="GroSeqPlus", hue="class", order=order, hue_order=hueOrder)
# add_stat_annotation(ax1, data=allDF, x="class", y="GroSeqPlus", hue="Condition", box_pairs=boxPairs, test='Mann-Whitney', loc='outside')
ax1.set(yscale="log")
plt.xticks(rotation=45)
plt.title("DMSO - GRO Plus Strand\nAnnotation Gene at Plus Strand")
plt.legend(loc="upper left")


fig.add_subplot(gs[1])
ax1 = sns.boxplot(data=dmsDF[dmsDF["geneStrand"] == 2], x="annotation", y="GroSeqPlus", hue="class", order=order, hue_order=hueOrder)
# add_stat_annotation(ax1, data=allDF, x="class", y="GroSeqPlus", hue="Condition", box_pairs=boxPairs, test='Mann-Whitney', loc='outside')
ax1.set(yscale="log")
plt.xticks(rotation=45)
plt.title("DMSO - GRO Plus Strand\nAnnotation Gene at Minus Strand")
plt.legend(loc="upper left")



fig.add_subplot(gs[2])
ax2 = sns.boxplot(data=dmsDF[dmsDF["geneStrand"] == 1], x="annotation", y="GroSeqMinus", hue="class", order=order, hue_order=hueOrder)
# add_stat_annotation(ax2, data=allDF, x="class", y="GroSeqMinus", hue="Condition", box_pairs=boxPairs, test='Mann-Whitney', loc='outside')
ax2.set(yscale="log")
plt.xticks(rotation=45)
plt.title("DMSO - GRO Minus Strand\nAnnotation Gene at Plus Strand")
plt.legend(loc="upper left")


fig.add_subplot(gs[3])
ax2 = sns.boxplot(data=dmsDF[dmsDF["geneStrand"] == 2], x="annotation", y="GroSeqMinus", hue="class", order=order, hue_order=hueOrder)
# add_stat_annotation(ax2, data=allDF, x="class", y="GroSeqMinus", hue="Condition", box_pairs=boxPairs, test='Mann-Whitney', loc='outside')
ax2.set(yscale="log")
plt.xticks(rotation=45)
plt.title("DMSO - GRO Minus Strand\nAnnotation Gene at Minus Strand")
plt.legend(loc="upper left")



fig.add_subplot(gs[4])
ax1 = sns.boxplot(data=dhtDF[dhtDF["geneStrand"] == 1], x="annotation", y="GroSeqPlus", hue="class", order=order, hue_order=hueOrder)
# add_stat_annotation(ax1, data=allDF, x="class", y="GroSeqPlus", hue="Condition", box_pairs=boxPairs, test='Mann-Whitney', loc='outside')
ax1.set(yscale="log")
plt.xticks(rotation=45)
plt.title("DHT - GRO Plus Strand\nAnnotation Gene at Plus Strand")
plt.legend(loc="upper left")



fig.add_subplot(gs[5])
ax1 = sns.boxplot(data=dhtDF[dhtDF["geneStrand"] == 2], x="annotation", y="GroSeqPlus", hue="class", order=order, hue_order=hueOrder)
# add_stat_annotation(ax1, data=allDF, x="class", y="GroSeqPlus", hue="Condition", box_pairs=boxPairs, test='Mann-Whitney', loc='outside')
ax1.set(yscale="log")
plt.xticks(rotation=45)
plt.title("DHT - GRO Plus Strand\nAnnotation Gene at Minus Strand")
plt.legend(loc="upper left")



fig.add_subplot(gs[6])
ax2 = sns.boxplot(data=dhtDF[dhtDF["geneStrand"] == 1], x="annotation", y="GroSeqMinus", hue="class", order=order, hue_order=hueOrder)
# add_stat_annotation(ax2, data=allDF, x="class", y="GroSeqMinus", hue="Condition", box_pairs=boxPairs, test='Mann-Whitney', loc='outside')
ax2.set(yscale="log")
plt.xticks(rotation=45)
plt.title("DHT - GRO Minus Strand\nAnnotation Gene at Plus Strand")
plt.legend(loc="upper left")



fig.add_subplot(gs[7])
ax2 = sns.boxplot(data=dhtDF[dhtDF["geneStrand"] == 2], x="annotation", y="GroSeqMinus", hue="class", order=order, hue_order=hueOrder)
# add_stat_annotation(ax2, data=allDF, x="class", y="GroSeqMinus", hue="Condition", box_pairs=boxPairs, test='Mann-Whitney', loc='outside')
ax2.set(yscale="log")
plt.xticks(rotation=45)
plt.title("DHT - GRO Minus Strand\nAnnotation Gene at Minus Strand")
plt.legend(loc="upper left")



fig.savefig(f"{figureRoot}/Gro.Annotations.all.pdf")

plt.close("all")



G = pickle.load(open(f"{dataRoot}/tmpData/GraphsGData.p", "rb" ))

# allDF = {"dht": dhtDF, "dms": dmsDF}
# pickle.dump(allDF, open(f"{dataRoot}/tmpData/allGRO.DF.p", "wb"))

allDF = pickle.load(open(f"{dataRoot}/tmpData/allGRO.DF.p", "rb" ))

dhtDF = allDF["dht"]
dmsDF = allDF["dms"]


dhtDF = dhtDF.sort_values("Gene").reset_index(drop=True)
dmsDF = dmsDF.sort_values("Gene").reset_index(drop=True)
dhtDF["mean"] = (dhtDF["GroSeqMinus"] + dhtDF["GroSeqPlus"]) / 2
dmsDF["mean"] = (dmsDF["GroSeqMinus"] + dmsDF["GroSeqPlus"]) / 2

dhtDF["sum"] = (dhtDF["GroSeqMinus"] + dhtDF["GroSeqPlus"])
dmsDF["sum"] = (dmsDF["GroSeqMinus"] + dmsDF["GroSeqPlus"])


[x for x in G.nodes() if x.find("KLK") != -1]

import time

dhtGroG = pd.DataFrame()
dmsGroG = pd.DataFrame()
i = 0

for edge in G.edges():
    if i % 100 == 0:
        t = time.time()
    aNode = edge[0]
    bNode = edge[1]
    left = bi.bisect_left(dhtDF["Gene"] , aNode)
    right = bi.bisect_right(dhtDF["Gene"] , aNode+"z")
    a = dhtDF[left:right]
    left = bi.bisect_left(dhtDF["Gene"] , bNode)
    right = bi.bisect_right(dhtDF["Gene"] , bNode+"z")
    b = dhtDF[left:right]
    if i % 100 == 0:
        print(f"first bisect: {time.time() - t}")
        t = time.time()
    if len(a) >= 1 and len(b) >= 1:
        a = a.rename(columns={"Gene": "aNode", "GroSeqMinus": "GroSeqMinusA", "GroSeqPlus": "GroSeqPlusA", "class": "classA", "annotation": "annotationA", "geneStrand": "geneStrandA", "Condition": "ConditionA"}).reset_index(drop=True)
        a = a.loc[a[["mean"]].idxmax(), :].reset_index(drop=True)
        a["nodeClassA"] = G.nodes[aNode]["nodeClass"]
        b = b.rename(columns={"Gene": "bNode", "GroSeqMinus": "GroSeqMinusB", "GroSeqPlus": "GroSeqPlusB", "class": "classB", "annotation": "annotationB", "geneStrand": "geneStrandB", "Condition": "ConditionB"}).reset_index(drop=True)
        b = b.loc[b[["mean"]].idxmax(), :].reset_index(drop=True)
        b["nodeClassB"] = G.nodes[bNode]["nodeClass"]
        ab = pd.concat([a, b], axis=1)
        dhtGroG = pd.concat([dhtGroG, ab], axis=0)
    if i % 100 == 0:
        print(f"first concat: {time.time() - t}")
        t = time.time()
    left = bi.bisect_left(dmsDF["Gene"] , aNode)
    right = bi.bisect_right(dmsDF["Gene"] , aNode+"z")
    a = dmsDF[left:right]
    left = bi.bisect_left(dmsDF["Gene"] , bNode)
    right = bi.bisect_right(dmsDF["Gene"] , bNode+"z")
    b = dmsDF[left:right]
    if i % 100 == 0:
        print(f"second bisect: {time.time() - t}")
        t = time.time()
    if len(a) >= 1 and len(b) >= 1:
        a = a.rename(columns={"Gene": "aNode", "GroSeqMinus": "GroSeqMinusA", "GroSeqPlus": "GroSeqPlusA", "class": "classA", "annotation": "annotationA", "geneStrand": "geneStrandA", "Condition": "ConditionA"}).reset_index(drop=True)
        a = a.loc[a[["mean"]].idxmax(), :].reset_index(drop=True)
        a["nodeClassA"] = G.nodes[aNode]["nodeClass"]
        b = b.rename(columns={"Gene": "bNode", "GroSeqMinus": "GroSeqMinusB", "GroSeqPlus": "GroSeqPlusB", "class": "classB", "annotation": "annotationB", "geneStrand": "geneStrandB", "Condition": "ConditionB"}).reset_index(drop=True)
        b = b.loc[b[["mean"]].idxmax(), :].reset_index(drop=True)
        b["nodeClassB"] = G.nodes[bNode]["nodeClass"]
        ab = pd.concat([a, b], axis=1)
        dmsGroG = pd.concat([dmsGroG, ab], axis=0)
    if i % 100 == 0:
        print(f"second concat: {time.time() - t}")
        t = time.time()
    if i % 100 == 0:
        print(i)
    i += 1


dmsGroG = dmsGroG.reset_index(drop=True)


for r in dmsGroG:
    aMinus = dmsGroG.loc[r, "GroSeqMinusA"]
    bMinus = dmsGroG.loc[r, "GroSeqMinusB"]
    aPlus = dmsGroG.loc[r, "GroSeqPlusA"]
    bPlus = dmsGroG.loc[r, "GroSeqPlusB"]
    aClass = dmsGroG.loc[r, "nodeClassA"]
    if dmsGroG.loc[r, "nodeClassA"] == "uPP":
        {}






allGroG = pd.concat([dmsGroG, dhtGroG])



fig = plt.figure(figsize=[15, 15])
gs = gridspec.GridSpec(ncols=2, nrows=2)
plt.subplots_adjust(bottom = 0.2, hspace=0.5)

fig.add_subplot(gs[0])
ax = sns.boxplot(y="GroSeqMinusA", x="ConditionB", hue="classA", data=allGroG, hue_order=hueOrder)
ax.set(yscale="log")
fig.add_subplot(gs[1])
ax = sns.boxplot(y="GroSeqMinusB", x="ConditionB", hue="classB", data=allGroG, hue_order=hueOrder)
ax.set(yscale="log")
fig.add_subplot(gs[2])
ax = sns.boxplot(y="GroSeqPlusA", x="ConditionB", hue="classA", data=allGroG, hue_order=hueOrder)
ax.set(yscale="log")
fig.add_subplot(gs[3])
ax = sns.boxplot(y="GroSeqPlusB", x="ConditionB", hue="classB", data=allGroG, hue_order=hueOrder)
ax.set(yscale="log")

fig.savefig(f"{figureRoot}/Gro.Edges.G.pdf")


import time

dhtGroN = pd.DataFrame()
dmsGroN = pd.DataFrame()
i = 0

for node in G.nodes():
    if i % 100 == 0:
        t = time.time()
    left = bi.bisect_left(dhtDF["Gene"] , node)
    right = bi.bisect_right(dhtDF["Gene"] , node+"z")
    a = dhtDF[left:right]
    a["degree"] = G.degree[node]
    a["nodClass"] = G.nodes[node]["nodeClass"]
    dhtGroN = pd.concat([dhtGroN, a], axis=0)
    if i % 100 == 0:
        print(f"first concat: {time.time() - t}")
        t = time.time()
    if i % 100 == 0:
        t = time.time()
    left = bi.bisect_left(dmsDF["Gene"] , node)
    right = bi.bisect_right(dmsDF["Gene"] , node+"z")
    a = dmsDF[left:right]
    a["degree"] = G.degree[node]
    a["nodClass"] = G.nodes[node]["nodeClass"]
    dmsGroN = pd.concat([dmsGroN, a], axis=0)
    if i % 100 == 0:
        print(f"second concat: {time.time() - t}")
        t = time.time()
        print(i)
    i += 1


nodeClass = ["con", "ind", "non"]

for node in nodeClass:
    scipy.stats.spearmanr(
        dhtGroN.loc[dhtGroN["nodClass"] == node, "sum"],
        dhtGroN.loc[dhtGroN["nodClass"] == node, "degree"]
        )


for node1 in nodeClass:
    print(node1)
    scipy.stats.ranksums(dhtGroN.loc[dhtGroN["nodClass"] == node1, "mean"], dmsGroN.loc[dmsGroN["nodClass"] == node1, "mean"])
    for node2 in nodeClass:
        print(node1, node2)
        scipy.stats.ranksums(dhtGroN.loc[dhtGroN["nodClass"] == node1, "mean"], dhtGroN.loc[dhtGroN["nodClass"] == node2, "mean"])
        scipy.stats.ranksums(dmsGroN.loc[dmsGroN["nodClass"] == node1, "mean"], dmsGroN.loc[dmsGroN["nodClass"] == node2, "mean"])





hueOrder=["con", "ind", "non", "oth", "tsP", "tsM"]
order = ["3' UTR", "5' UTR", "Distal Intergenic", "Downstream ",  "Exon ", "Intron ", "Promoter "]

dmsGroN = dmsGroN.sort_values("degree").reset_index(drop=True)
dhtGroN = dhtGroN.sort_values("degree").reset_index(drop=True)


nodeClasses = ["con", "ind", "non", "oth", "tsP", "tsM", "uPP", "dPP", "uMP", "dMP"]
for node in nodeClasses:
    df = dmsGroN[dmsGroN["nodClass"] == node]
    r, p = scipy.stats.spearmanr(df["degree"], df["GroSeqMinus"])
    print(f"Gro Minus {node} class has a spearman correlation coefficient {round(r, 3)} with {round(p, 3)} p-value from DMSO.")
    df = dhtGroN[dhtGroN["nodClass"] == node]
    r, p = scipy.stats.spearmanr(df["degree"], df["GroSeqMinus"])
    print(f"Gro Minus {node} class has a spearman correlation coefficient {round(r, 3)} with {round(p, 3)} p-value from DHT.")
    df = dmsGroN[dmsGroN["nodClass"] == node]
    r, p = scipy.stats.spearmanr(df["degree"], df["GroSeqPlus"])
    print(f"Gro Plus {node} class has a spearman correlation coefficient {round(r, 3)} with {round(p, 3)} p-value from DMSO.")
    df = dhtGroN[dhtGroN["nodClass"] == node]
    r, p = scipy.stats.spearmanr(df["degree"], df["GroSeqPlus"])
    print(f"Gro Plus {node} class has a spearman correlation coefficient {round(r, 3)} with {round(p, 3)} p-value from DHT.")


fig = plt.figure(figsize=[15, 15])

g = sns.scatterplot(x="degree", y="GroSeqPlus", hue="class", data=dmsGroN)
# g = sns.lmplot(x="degree", y="GroSeqMinus", hue="class", col="geneStrand", data=dmsGroN)
fig.savefig(f"{figureRoot}/Gro.Degree.G.lm.pdf")


fig = plt.figure(figsize=[15, 15])
gs = gridspec.GridSpec(ncols=2, nrows=2)
plt.subplots_adjust(bottom = 0.2, hspace=0.5)





fig.add_subplot(gs[0])
ax = sns.lmplot(x="degree", y="GroSeqPlus", hue="class", data=dmsGroN)
ax.set(yscale="log")
plt.title("DMSO - Gro Plus")

fig.add_subplot(gs[1])
ax = sns.lmplot(x="degree", y="GroSeqMinus", hue="class", data=dmsGroN, hue_order=hueOrder)
ax.set(yscale="log")
plt.title("DMSO - Gro Minus")

fig.add_subplot(gs[2])
ax = sns.lmplot(x="degree", y="GroSeqPlus", hue="class", data=dhtGroN, hue_order=hueOrder)
ax.set(yscale="log")
plt.title("DHT - Gro Plus")

fig.add_subplot(gs[3])
ax = sns.lmplot(x="degree", y="GroSeqMinus", hue="class", data=dhtGroN, hue_order=hueOrder)
ax.set(yscale="log")
plt.title("DHT - Gro Minus")

fig.savefig(f"{figureRoot}/Gro.Degree.G.lm.pdf")


df = dhtGroG.groupby(["annotationA", "annotationB"]).size().reset_index()
df = dhtGroG.groupby(["classA", "classB"]).size().reset_index()

df2 = df.pivot(index="annotationA", columns="annotationB", values=0).fillna(0)
df2 = df.pivot(index="classA", columns="classB", values=0).fillna(0)


from Functions.PlotFunctions import *

hc = ["#BACDCC", "#5DBACD", "#3688B2", "#1A6588", "#1A3249"]

th = [0, 0.15, 0.5, 0.85, 1]
cdict = NonLinCdict(th, hc)
cm = colors.LinearSegmentedColormap('test', cdict)




sizes = [318,134,7887,99,734,8387,20128]

matIF = []
for i in range(6):
    for j in range(6):
        if i == j:
            V = sizes[i]
            Emax = (V*(V+1))/2
        else:
            U = sizes[i]
            V = sizes[j]
            Emax = U*V
        IF = (df2.iloc[i,j] + df2.iloc[j,i]) /Emax*2
        matIF += [IF]

matIF = np.array(matIF).reshape(6,6)


allGro = pd.concat([dhtGroG, dmsGroG])
fig = plt.figure(figsize=[6, 6])
gs = gridspec.GridSpec(ncols=1, nrows=1)
plt.subplots_adjust(bottom = 0.2, hspace=0.5)

sns.heatmap(matIF, cmap=cm, xticklabels=order, yticklabels=order)

fig.savefig(f"{figureRoot}/Annotation.Interaction.G.heat.pdf")


dhtGroG
