

from Functions.Packages import *
from Functions.PlotFunctions import *

G = pickle.load(open(f"{dataRoot}/tmpData/GraphsGData.p", "rb" ))


veh = pd.read_csv(f"{dataRoot}/LNCaP.PRO.Cap/GSM1543793_LNCaP_PRO-cap_1h_vehicle_hg19_peaks.bed", sep="\t", names=["chr.veh", "start.veh", "end.veh", "score.veh", "signal.veh", "strand.veh"])
dht = pd.read_csv(f"{dataRoot}/LNCaP.PRO.Cap/GSM1543794_LNCaP_PRO-cap_1h_dht_hg19_peaks.bed", sep="\t", names=["chr.dht", "start.dht", "end.dht", "score.dht", "signal.dht", "strand.dht"])

conAnn = pd.read_csv("/home/ualtintas/ARBSs/annotations/con.Annotations.txt", sep="\t")
conAnn = conAnn.rename(columns={"PeakID (cmd=annotatePeaks.pl regions/cons-arbs.bed hg19)": "PeakID"})
indAnn = pd.read_csv("/home/ualtintas/ARBSs/annotations/ind.Annotations.txt", sep="\t")
indAnn = indAnn.rename(columns={"PeakID (cmd=annotatePeaks.pl regions/ind-arbs.bed hg19)": "PeakID"})
nonAnn = pd.read_csv("/home/ualtintas/ARBSs/annotations/non.Annotations.txt", sep="\t")
nonAnn = nonAnn.rename(columns={"PeakID (cmd=annotatePeaks.pl regions/Non-Active-ARBS.bed hg19)": "PeakID"})




conAnn = pd.read_csv(f"{dataRoot}/Regions/con.AnnFile.cvs", sep="\t")
indAnn = pd.read_csv(f"{dataRoot}/Regions/ind.AnnFile.csv", sep="\t")
nonAnn = pd.read_csv(f"{dataRoot}/Regions/non.AnnFile.csv", sep="\t")
tsPAnn = pd.read_csv("/home/ualtintas/genomeAnnotations/Regions/TSS.hg19.+.AnnFile.csv", sep="\t")
tsMAnn = pd.read_csv("/home/ualtintas/genomeAnnotations/Regions/TSS.hg19.-.AnnFile.csv", sep="\t")

conAnn = conAnn.rename(columns={"seqnames": "Chr", "start": "Start", "end": "End"})
indAnn = indAnn.rename(columns={"seqnames": "Chr", "start": "Start", "end": "End"})
nonAnn = nonAnn.rename(columns={"seqnames": "Chr", "start": "Start", "end": "End"})
tsPAnn = tsPAnn.rename(columns={"seqnames": "Chr", "start": "Start", "end": "End"})
tsMAnn = tsMAnn.rename(columns={"seqnames": "Chr", "start": "Start", "end": "End"})


veh = veh.sort_values("start.veh").reset_index(drop=True)
dht = dht.sort_values("start.dht").reset_index(drop=True)


# capG = capG.sort_values("SYMBOL").reset_index(drop=True)

logFC = pd.read_csv(f"{dataRoot}/DEG/GSE64529_diffexpr-results.csv")
logFC = logFC[["Gene", "log2FoldChange", "padj"]]
logFC = logFC.sort_values("Gene").reset_index(drop=True)



conAnn["SignalVEH.+"] = 0
conAnn["SignalDHT.+"] = 0


conAnn["SignalVEH.-"] = 0
conAnn["SignalDHT.-"] = 0

conAnn["SignalVEH"] = 0
conAnn["SignalDHT"] = 0

conAnn["NumberofCapsVEH.+"] = 0
conAnn["NumberofCapsDHT.+"] = 0


conAnn["NumberofCapsVEH.-"] = 0
conAnn["NumberofCapsDHT.-"] = 0

conAnn["NumberofCapsVEH"] = 0
conAnn["NumberofCapsDHT"] = 0



for i in range(len(conAnn)):
    chr = conAnn.loc[i, "Chr"]
    start = conAnn.loc[i, "Start"]
    end = conAnn.loc[i, "End"]
    sym = conAnn.loc[i, "SYMBOL"]
    leftV = bi.bisect_left(veh["start.veh"] , start)
    rightV = bi.bisect_left(veh["end.veh"] , end)
    leftD = bi.bisect_left(dht["start.dht"] , start)
    rightD = bi.bisect_left(dht["end.dht"] , end)
    tmpV = veh[leftV:rightV]
    tmpV = tmpV[tmpV["chr.veh"] == chr]
    if tmpV.shape[0] > 0:
        conAnn.loc[i, "SignalVEH.+"] = max(max(list(tmpV[tmpV["strand.veh"] == "+"]["signal.veh"]), [0]))
        conAnn.loc[i, "SignalVEH.-"] = max(max(list(tmpV[tmpV["strand.veh"] == "-"]["signal.veh"]), [0]))
        conAnn.loc[i, "NumberofCapsVEH.+"] = tmpV[tmpV["strand.veh"] == "+"].shape[0]
        conAnn.loc[i, "NumberofCapsVEH.-"] = tmpV[tmpV["strand.veh"] == "-"].shape[0]
        conAnn.loc[i,"SignalVEH"] = tmpV["signal.veh"].max()
        conAnn.loc[i, "NumberofCapsVEH"] = tmpV.shape[0]
    tmpD = dht[leftD:rightD]
    tmpD = tmpD[tmpD["chr.dht"] == chr]
    if tmpD.shape[0] > 0:
        conAnn.loc[i, "SignalDHT.+"] = max(max(list(tmpD[tmpD["strand.dht"] == "+"]["signal.dht"]), [0]))
        conAnn.loc[i, "SignalDHT.-"] = max(max(list(tmpD[tmpD["strand.dht"] == "-"]["signal.dht"]), [0]))
        conAnn.loc[i, "NumberofCapsDHT.+"] = tmpD[tmpD["strand.dht"] == "+"].shape[0]
        conAnn.loc[i, "NumberofCapsDHT.-"] = tmpD[tmpD["strand.dht"] == "-"].shape[0]
        conAnn.loc[i,"SignalDHT"] = tmpD["signal.dht"].max()
        conAnn.loc[i, "NumberofCaplDHT"] = tmpD.shape[0]





indAnn["SignalVEH.+"] = 0
indAnn["SignalDHT.+"] = 0

indAnn["SignalVEH.-"] = 0
indAnn["SignalDHT.-"] = 0

indAnn["SignalVEH"] = 0
indAnn["SignalDHT"] = 0

indAnn["NumberofCapsVEH.+"] = 0
indAnn["NumberofCapsDHT.+"] = 0

indAnn["NumberofCapsVEH.-"] = 0
indAnn["NumberofCapsDHT.-"] = 0

indAnn["NumberofCapsVEH"] = 0
indAnn["NumberofCapsDHT"] = 0



for i in range(len(indAnn)):
    chr = indAnn.loc[i, "Chr"]
    start = indAnn.loc[i, "Start"]
    end = indAnn.loc[i, "End"]
    sym = indAnn.loc[i, "SYMBOL"]
    leftV = bi.bisect_left(veh["start.veh"] , start)
    rightV = bi.bisect_left(veh["end.veh"] , end)
    leftD = bi.bisect_left(dht["start.dht"] , start)
    rightD = bi.bisect_left(dht["end.dht"] , end)
    tmpV = veh[leftV:rightV]
    tmpV = tmpV[tmpV["chr.veh"] == chr]
    if tmpV.shape[0] > 0:
        indAnn.loc[i, "SignalVEH.+"] = max(max(list(tmpV[tmpV["strand.veh"] == "+"]["signal.veh"]), [0]))
        indAnn.loc[i, "SignalVEH.-"] = max(max(list(tmpV[tmpV["strand.veh"] == "-"]["signal.veh"]), [0]))
        indAnn.loc[i, "NumberofCapsVEH.+"] = tmpV[tmpV["strand.veh"] == "+"].shape[0]
        indAnn.loc[i, "NumberofCapsVEH.-"] = tmpV[tmpV["strand.veh"] == "-"].shape[0]
        indAnn.loc[i,"SignalVEH"] = tmpV["signal.veh"].max()
        indAnn.loc[i, "NumberofCapsVEH"] = tmpV.shape[0]
    tmpD = dht[leftD:rightD]
    tmpD = tmpD[tmpD["chr.dht"] == chr]
    if tmpD.shape[0] > 0:
        indAnn.loc[i, "SignalDHT.+"] = max(max(list(tmpD[tmpD["strand.dht"] == "+"]["signal.dht"]), [0]))
        indAnn.loc[i, "SignalDHT.-"] = max(max(list(tmpD[tmpD["strand.dht"] == "-"]["signal.dht"]), [0]))
        indAnn.loc[i, "NumberofCapsDHT.+"] = tmpD[tmpD["strand.dht"] == "+"].shape[0]
        indAnn.loc[i, "NumberofCapsDHT.-"] = tmpD[tmpD["strand.dht"] == "-"].shape[0]
        indAnn.loc[i,"SignalDHT"] = tmpD["signal.dht"].max()
        indAnn.loc[i, "NumberofCaplDHT"] = tmpD.shape[0]




nonAnn["SignalVEH.+"] = 0
nonAnn["SignalDHT.+"] = 0

nonAnn["SignalVEH.-"] = 0
nonAnn["SignalDHT.-"] = 0

nonAnn["SignalVEH"] = 0
nonAnn["SignalDHT"] = 0

nonAnn["NumberofCapsVEH.+"] = 0
nonAnn["NumberofCapsDHT.+"] = 0

nonAnn["NumberofCapsVEH.-"] = 0
nonAnn["NumberofCapsDHT.-"] = 0

nonAnn["NumberofCapsVEH"] = 0
nonAnn["NumberofCapsDHT"] = 0




for i in range(len(nonAnn)):
    chr = nonAnn.loc[i, "Chr"]
    start = nonAnn.loc[i, "Start"]
    end = nonAnn.loc[i, "End"]
    leftV = bi.bisect_left(veh["start.veh"] , start)
    rightV = bi.bisect_left(veh["end.veh"] , end)
    leftD = bi.bisect_left(dht["start.dht"] , start)
    rightD = bi.bisect_left(dht["end.dht"] , end)
    tmpV = veh[leftV:rightV]
    tmpV = tmpV[tmpV["chr.veh"] == chr]
    if tmpV.shape[0] > 0:
        nonAnn.loc[i, "SignalVEH.+"] = max(max(list(tmpV[tmpV["strand.veh"] == "+"]["signal.veh"]), [0]))
        nonAnn.loc[i, "SignalVEH.-"] = max(max(list(tmpV[tmpV["strand.veh"] == "-"]["signal.veh"]), [0]))
        nonAnn.loc[i, "NumberofCapsVEH.+"] = tmpV[tmpV["strand.veh"] == "+"].shape[0]
        nonAnn.loc[i, "NumberofCapsVEH.-"] = tmpV[tmpV["strand.veh"] == "-"].shape[0]
        nonAnn.loc[i,"SignalVEH"] = tmpV["signal.veh"].max()
        nonAnn.loc[i, "NumberofCapsVEH"] = tmpV.shape[0]
    tmpD = dht[leftD:rightD]
    tmpD = tmpD[tmpD["chr.dht"] == chr]
    if tmpD.shape[0] > 0:
        nonAnn.loc[i, "SignalDHT.+"] = max(max(list(tmpD[tmpD["strand.dht"] == "+"]["signal.dht"]), [0]))
        nonAnn.loc[i, "SignalDHT.-"] = max(max(list(tmpD[tmpD["strand.dht"] == "-"]["signal.dht"]), [0]))
        nonAnn.loc[i, "NumberofCapsDHT.+"] = tmpD[tmpD["strand.dht"] == "+"].shape[0]
        nonAnn.loc[i, "NumberofCapsDHT.-"] = tmpD[tmpD["strand.dht"] == "-"].shape[0]
        nonAnn.loc[i,"SignalDHT"] = tmpD["signal.dht"].max()
        nonAnn.loc[i, "NumberofCaplDHT"] = tmpD.shape[0]



tsPAnn["SignalVEH.+"] = 0
tsPAnn["SignalDHT.+"] = 0

tsPAnn["SignalVEH.-"] = 0
tsPAnn["SignalDHT.-"] = 0

tsPAnn["SignalVEH"] = 0
tsPAnn["SignalDHT"] = 0

tsPAnn["NumberofCapsVEH.+"] = 0
tsPAnn["NumberofCapsDHT.+"] = 0

tsPAnn["NumberofCapsVEH.-"] = 0
tsPAnn["NumberofCapsDHT.-"] = 0

tsPAnn["NumberofCapsVEH"] = 0
tsPAnn["NumberofCapsDHT"] = 0




for i in range(len(tsPAnn)):
    chr = tsPAnn.loc[i, "Chr"]
    start = tsPAnn.loc[i, "Start"]
    end = tsPAnn.loc[i, "End"]
    leftV = bi.bisect_left(veh["start.veh"] , start)
    rightV = bi.bisect_left(veh["end.veh"] , end)
    leftD = bi.bisect_left(dht["start.dht"] , start)
    rightD = bi.bisect_left(dht["end.dht"] , end)
    tmpV = veh[leftV:rightV]
    tmpV = tmpV[tmpV["chr.veh"] == chr]
    if tmpV.shape[0] > 0:
        tsPAnn.loc[i, "SignalVEH.+"] = max(max(list(tmpV[tmpV["strand.veh"] == "+"]["signal.veh"]), [0]))
        tsPAnn.loc[i, "SignalVEH.-"] = max(max(list(tmpV[tmpV["strand.veh"] == "-"]["signal.veh"]), [0]))
        tsPAnn.loc[i, "NumberofCapsVEH.+"] = tmpV[tmpV["strand.veh"] == "+"].shape[0]
        tsPAnn.loc[i, "NumberofCapsVEH.-"] = tmpV[tmpV["strand.veh"] == "-"].shape[0]
        tsPAnn.loc[i,"SignalVEH"] = tmpV["signal.veh"].max()
        tsPAnn.loc[i, "NumberofCapsVEH"] = tmpV.shape[0]
    tmpD = dht[leftD:rightD]
    tmpD = tmpD[tmpD["chr.dht"] == chr]
    if tmpD.shape[0] > 0:
        tsPAnn.loc[i, "SignalDHT.+"] = max(max(list(tmpD[tmpD["strand.dht"] == "+"]["signal.dht"]), [0]))
        tsPAnn.loc[i, "SignalDHT.-"] = max(max(list(tmpD[tmpD["strand.dht"] == "-"]["signal.dht"]), [0]))
        tsPAnn.loc[i, "NumberofCapsDHT.+"] = tmpD[tmpD["strand.dht"] == "+"].shape[0]
        tsPAnn.loc[i, "NumberofCapsDHT.-"] = tmpD[tmpD["strand.dht"] == "-"].shape[0]
        tsPAnn.loc[i,"SignalDHT"] = tmpD["signal.dht"].max()
        tsPAnn.loc[i, "NumberofCaplDHT"] = tmpD.shape[0]


tsMAnn["SignalVEH.+"] = 0
tsMAnn["SignalDHT.+"] = 0

tsMAnn["SignalVEH.-"] = 0
tsMAnn["SignalDHT.-"] = 0

tsMAnn["SignalVEH"] = 0
tsMAnn["SignalDHT"] = 0

tsMAnn["NumberofCapsVEH.+"] = 0
tsMAnn["NumberofCapsDHT.+"] = 0

tsMAnn["NumberofCapsVEH.-"] = 0
tsMAnn["NumberofCapsDHT.-"] = 0

tsMAnn["NumberofCapsVEH"] = 0
tsMAnn["NumberofCapsDHT"] = 0




for i in range(len(tsMAnn)):
    chr = tsMAnn.loc[i, "Chr"]
    start = tsMAnn.loc[i, "Start"]
    end = tsMAnn.loc[i, "End"]
    leftV = bi.bisect_left(veh["start.veh"] , start)
    rightV = bi.bisect_left(veh["end.veh"] , end)
    leftD = bi.bisect_left(dht["start.dht"] , start)
    rightD = bi.bisect_left(dht["end.dht"] , end)
    tmpV = veh[leftV:rightV]
    tmpV = tmpV[tmpV["chr.veh"] == chr]
    if tmpV.shape[0] > 0:
        tsMAnn.loc[i, "SignalVEH.+"] = max(max(list(tmpV[tmpV["strand.veh"] == "+"]["signal.veh"]), [0]))
        tsMAnn.loc[i, "SignalVEH.-"] = max(max(list(tmpV[tmpV["strand.veh"] == "-"]["signal.veh"]), [0]))
        tsMAnn.loc[i, "NumberofCapsVEH.+"] = tmpV[tmpV["strand.veh"] == "+"].shape[0]
        tsMAnn.loc[i, "NumberofCapsVEH.-"] = tmpV[tmpV["strand.veh"] == "-"].shape[0]
        tsMAnn.loc[i,"SignalVEH"] = tmpV["signal.veh"].max()
        tsMAnn.loc[i, "NumberofCapsVEH"] = tmpV.shape[0]
    tmpD = dht[leftD:rightD]
    tmpD = tmpD[tmpD["chr.dht"] == chr]
    if tmpD.shape[0] > 0:
        tsMAnn.loc[i, "SignalDHT.+"] = max(max(list(tmpD[tmpD["strand.dht"] == "+"]["signal.dht"]), [0]))
        tsMAnn.loc[i, "SignalDHT.-"] = max(max(list(tmpD[tmpD["strand.dht"] == "-"]["signal.dht"]), [0]))
        tsMAnn.loc[i, "NumberofCapsDHT.+"] = tmpD[tmpD["strand.dht"] == "+"].shape[0]
        tsMAnn.loc[i, "NumberofCapsDHT.-"] = tmpD[tmpD["strand.dht"] == "-"].shape[0]
        tsMAnn.loc[i,"SignalDHT"] = tmpD["signal.dht"].max()
        tsMAnn.loc[i, "NumberofCaplDHT"] = tmpD.shape[0]






con = conAnn[["V4", "Chr", "Start", "End", "annotation", "geneStrand", "SignalVEH.-", "NumberofCapsVEH.-", "SignalVEH.+", "NumberofCapsVEH.+", "SignalVEH", "NumberofCapsVEH", "SignalDHT.-", "NumberofCapsDHT.-", "SignalDHT.+", "NumberofCapsDHT.+", "SignalDHT", "NumberofCapsDHT", "SYMBOL"]]
ind = indAnn[["V4", "Chr", "Start", "End", "annotation", "geneStrand", "SignalVEH.-", "NumberofCapsVEH.-", "SignalVEH.+", "NumberofCapsVEH.+", "SignalVEH", "NumberofCapsVEH", "SignalDHT.-", "NumberofCapsDHT.-", "SignalDHT.+", "NumberofCapsDHT.+", "SignalDHT", "NumberofCapsDHT", "SYMBOL"]]
non = nonAnn[["V4", "Chr", "Start", "End", "annotation", "geneStrand", "SignalVEH.-", "NumberofCapsVEH.-", "SignalVEH.+", "NumberofCapsVEH.+", "SignalVEH", "NumberofCapsVEH", "SignalDHT.-", "NumberofCapsDHT.-", "SignalDHT.+", "NumberofCapsDHT.+", "SignalDHT", "NumberofCapsDHT", "SYMBOL"]]
tsP = tsPAnn[["V4", "Chr", "Start", "End", "annotation", "geneStrand", "SignalVEH.-", "NumberofCapsVEH.-", "SignalVEH.+", "NumberofCapsVEH.+", "SignalVEH", "NumberofCapsVEH", "SignalDHT.-", "NumberofCapsDHT.-", "SignalDHT.+", "NumberofCapsDHT.+", "SignalDHT", "NumberofCapsDHT", "SYMBOL"]]
tsM = tsMAnn[["V4", "Chr", "Start", "End", "annotation", "geneStrand", "SignalVEH.-", "NumberofCapsVEH.-", "SignalVEH.+", "NumberofCapsVEH.+", "SignalVEH", "NumberofCapsVEH", "SignalDHT.-", "NumberofCapsDHT.-", "SignalDHT.+", "NumberofCapsDHT.+", "SignalDHT", "NumberofCapsDHT", "SYMBOL"]]


con["nodeClass"] = "con"
ind["nodeClass"] = "ind"
non["nodeClass"] = "non"
tsP["nodeClass"] = "tsP"
tsM["nodeClass"] = "tsM"


capG = pd.concat([con, ind, non, tsP, tsM]).reset_index(drop=True)

capG["annotation"] = capG["annotation"].str.split("(", expand=True)[0]

df = capG.groupby(["nodeClass", "annotation"]).size().reset_index()

df.pivot(index="annotation", columns="nodeClass", values=0)



capG["SignalVEH.-"] += 1
capG["SignalVEH.+"] += 1

capG["SignalDHT.-"] += 1
capG["SignalDHT.+"] += 1


capG["annotation2"] = "other"
capG.loc[(capG["annotation"].isin(["Promoter "])), "annotation2"] = "Promoter"


df = capG.groupby(["nodeClass", "annotation2"]).mean().reset_index()

df.pivot(index="annotation2", columns="nodeClass", values="SignalDHT.-")










nodeClasses = ["con", "ind", "non", "tsP", "tsM"]
for node in nodeClasses:
    df = capG[capG["nodeClass"] == node]
    proDF = df[df["annotation2"] == "Promoter"]
    othDF = df[df["annotation2"] == "other"]
    r, p = scipy.stats.ranksums(proDF["SignalVEH.-"], othDF["SignalVEH.-"])
    print(f"Pro Minus - VEH {node} class's annotated promoters vs others, Wilcoxon rank-sum test statistic:{round(r, 3)}, with pVal:{round(p, 3)}.")
    r, p = scipy.stats.ranksums(proDF["SignalVEH.+"], othDF["SignalVEH.+"])
    print(f"Pro Plus - VEH {node} class's annotated promoters vs others, Wilcoxon rank-sum test statistic:{round(r, 3)}, with pVal:{round(p, 3)}.")
    r, p = scipy.stats.ranksums(proDF["SignalDHT.-"], othDF["SignalDHT.-"])
    print(f"Pro Minus - DHT {node} class's annotated promoters vs others, Wilcoxon rank-sum test statistic:{round(r, 3)}, with pVal:{round(p, 3)}.")
    r, p = scipy.stats.ranksums(proDF["SignalDHT.+"], othDF["SignalDHT.+"])
    print(f"Pro Plus - DHT {node} class's annotated promoters vs others, Wilcoxon rank-sum test statistic:{round(r, 3)}, with pVal:{round(p, 3)}.")



pickle.dump(capG, open("capG.p", "wb"))


order=["con", "ind", "non", "tsP", "tsM"]
hueOrder = ["Promoter", "other"]

fig = plt.figure(figsize=[15, 15])
gs = gridspec.GridSpec(ncols=2, nrows=2)
plt.subplots_adjust(bottom = 0.2, hspace=0.5)


fig.add_subplot(gs[0])
ax1 = sns.violinplot(data=capG, hue="annotation2", y="SignalVEH.+", x="nodeClass", order=order, hue_order=hueOrder, split=False)
# add_stat_annotation(ax1, data=allDF, x="class", y="GroSeqPlus", hue="Condition", box_pairs=boxPairs, test='Mann-Whitney', loc='outside')
# ax1.set(yscale="log")
plt.xticks(rotation=45)
plt.title("VEH - Pro Cap Plus Signal")
plt.legend(loc="upper left")


fig.add_subplot(gs[1])
ax1 = sns.violinplot(data=capG, hue="annotation2", y="SignalVEH.-", x="nodeClass", order=order, hue_order=hueOrder, split=False)
# add_stat_annotation(ax1, data=allDF, x="class", y="GroSeqPlus", hue="Condition", box_pairs=boxPairs, test='Mann-Whitney', loc='outside')
# ax1.set(yscale="log")
plt.xticks(rotation=45)
plt.title("VEH - Pro Cap Minus Signal")
plt.legend(loc="upper left")



fig.add_subplot(gs[2])
ax2 = sns.violinplot(data=capG, hue="annotation2", y="SignalDHT.+", x="nodeClass", order=order, hue_order=hueOrder, split=False)
# add_stat_annotation(ax2, data=allDF, x="class", y="GroSeqMinus", hue="Condition", box_pairs=boxPairs, test='Mann-Whitney', loc='outside')
# ax2.set(yscale="log")
plt.xticks(rotation=45)
plt.title("DHT - Pro Cap Plus Signal")
plt.legend(loc="upper left")


fig.add_subplot(gs[3])
ax2 = sns.violinplot(data=capG, hue="annotation2", y="SignalDHT.-", x="nodeClass", order=order, hue_order=hueOrder, split=False)
# add_stat_annotation(ax2, data=allDF, x="class", y="GroSeqMinus", hue="Condition", box_pairs=boxPairs, test='Mann-Whitney', loc='outside')
# ax2.set(yscale="log")
plt.xticks(rotation=45)
plt.title("DHT - Pro Cap Minus Signal")
plt.legend(loc="upper left")




fig.savefig(f"{figureRoot}/Pro.AnnotationsVsOthers.all.pdf")




# RainCloudPlot


import matplotlib.collections as clt


import ptitprince as pt

dx = "nodeClass"



dx = "nodeClass"; dhue = "annotation2"; ort = "v"; pal = "Set2"; sigma = .2



fig = plt.figure(figsize=[15, 15])
gs = gridspec.GridSpec(ncols=2, nrows=2)
plt.subplots_adjust(bottom = 0.2, hspace=0.5)


ax1 = fig.add_subplot(gs[0])

dy = "SignalVEH.+"
ax1 = pt.RainCloud(x = dx, y = dy, hue = dhue, data = capG, palette = pal, bw = sigma, width_viol = .7,
                ax = ax1, orient = ort , alpha = .65, dodge = True, move = .2)
ax1.get_legend().remove()
ax1.set(yscale="log")

ax2 = fig.add_subplot(gs[1])

dy = "SignalVEH.-"
ax2 = pt.RainCloud(x = dx, y = dy, hue = dhue, data = capG, palette = pal, bw = sigma, width_viol = .7,
                ax = ax2, orient = ort , alpha = .65, dodge = True, move = .2)

ax2.set(yscale="log")


ax3 = fig.add_subplot(gs[2])

dy = "SignalDHT.+";
ax3 = pt.RainCloud(x = dx, y = dy, hue = dhue, data = capG, palette = pal, bw = sigma, width_viol = .7,
                ax = ax3, orient = ort , alpha = .65, dodge = True, move = .2)
ax3.get_legend().remove()
ax3.set(yscale="log")

ax4 = fig.add_subplot(gs[3])

dy = "SignalDHT.-";
ax4 = pt.RainCloud(x = dx, y = dy, hue = dhue, data = capG, palette = pal, bw = sigma, width_viol = .7,
                ax = ax4, orient = ort , alpha = .65, dodge = True, move = .2)
ax4.get_legend().remove()
ax4.set(yscale="log")

fig.savefig(f"{figureRoot}/Pro.AnnotationsVsOthers.all.pdf")












ax=pt.half_violinplot( x = dx, y = "SignalVEH.+", data = df, palette = pal, bw = .2, cut = 0.,
                      scale = "area", width = .6, inner = None, orient = ort)

ax=sns.stripplot(x = dx, y = dy, data = df, palette = pal, edgecolor = "white",
                 size = 3, jitter = 1, zorder = 0, orient = ort)

ax=sns.boxplot( x = dx, y = dy, data = df, color = "black", width = .15, zorder = 10,\
            showcaps = True, boxprops = {'facecolor':'none', "zorder":10},\
            showfliers=True, whiskerprops = {'linewidth':2, "zorder":10},\
               saturation = 1, orient = ort)










fig.add_subplot(gs[4])
ax1 = sns.violinplot(data=capG[capG["geneStrand"] == 1], x="annotation", y="SignalDHT.+", hue="nodeClass", order=order, hue_order=hueOrder, inner="stick", scale_hue=False, split=True)
# add_stat_annotation(ax1, data=allDF, x="class", y="GroSeqPlus", hue="Condition", box_pairs=boxPairs, test='Mann-Whitney', loc='outside')
ax1.set(yscale="log")
plt.xticks(rotation=45)
plt.title("DHT - Pro Cap Plus Signal\nAnnotation Gene at Plus Strand")
plt.legend(loc="upper left")



fig.add_subplot(gs[5])
ax1 = sns.boxplot(data=capG[capG["geneStrand"] == 2], x="annotation", y="SignalDHT.+", hue="nodeClass", order=order, hue_order=hueOrder, inner="stick", scale_hue=False, split=True)
# add_stat_annotation(ax1, data=allDF, x="class", y="GroSeqPlus", hue="Condition", box_pairs=boxPairs, test='Mann-Whitney', loc='outside')
ax1.set(yscale="log")
plt.xticks(rotation=45)
plt.title("DHT - Pro Cap Plus Signal\nAnnotation Gene at Minus Strand")
plt.legend(loc="upper left")



fig.add_subplot(gs[6])
ax2 = sns.boxplot(data=capG[capG["geneStrand"] == 1], x="annotation", y="SignalDHT.-", hue="nodeClass", order=order, hue_order=hueOrder, inner="stick", scale_hue=False, split=True)
# add_stat_annotation(ax2, data=allDF, x="class", y="GroSeqMinus", hue="Condition", box_pairs=boxPairs, test='Mann-Whitney', loc='outside')
ax2.set(yscale="log")
plt.xticks(rotation=45)
plt.title("DHT - Pro Cap Minus Signal\nAnnotation Gene at Plus Strand")
plt.legend(loc="upper left")



fig.add_subplot(gs[7])
ax2 = sns.boxplot(data=capG[capG["geneStrand"] == 2], x="annotation", y="SignalDHT.-", hue="nodeClass", order=order, hue_order=hueOrder, inner="stick", scale_hue=False, split=True)
# add_stat_annotation(ax2, data=allDF, x="class", y="GroSeqMinus", hue="Condition", box_pairs=boxPairs, test='Mann-Whitney', loc='outside')
ax2.set(yscale="log")
plt.xticks(rotation=45)
plt.title("DHT - Pro Cap Minus Signal\nAnnotation Gene at Minus Strand")
plt.legend(loc="upper left")



plt.close("all")


capSig = capG[["nodeClass", "SignalDHT.+", "SignalDHT.-", "SignalVEH.+", "SignalVEH.-", "V4"]]


dhtCap = capSig[["nodeClass", "SignalDHT.+", "SignalDHT.-", "V4"]]
dhtCap.loc[:, "Condition"] = "DHT"
dhtCap = dhtCap.rename(columns={"SignalDHT.+": "Signal.+", "SignalDHT.-": "Signal.-"})

vehCap = capSig[["nodeClass", "SignalVEH.+", "SignalVEH.-", "V4"]]
vehCap.loc[:,"Condition"] = "VEH"
vehCap = vehCap.rename(columns={"SignalVEH.+": "Signal.+", "SignalVEH.-": "Signal.-"})

allCap = pd.concat([dhtCap,vehCap])


allCap[(allCap["Signal.-"] > 600) & (allCap["Signal.+"]< 600)]



hueOrder=["con", "ind", "non", "oth", "tsP", "tsM"]
order = ["3' UTR", "5' UTR", "Distal Intergenic", "Downstream ",  "Exon ", "Intron ", "Promoter "]
colorPalette = {"VEH": "#FAB550", "DHT": "#5F9BD2"}


fig = plt.figure(figsize=[15,6])
plt.subplots_adjust(bottom = 0.2)
gs = gridspec.GridSpec(ncols=2, nrows=1)
fig.add_subplot(gs[0])
ax1 = sns.barplot(hue="Condition", y="Signal.-", x="nodeClass", data=allCap, order=hueOrder, hue_order=["VEH", "DHT"], palette=colorPalette)
plt.xticks(rotation=45)
plt.legend()

fig.add_subplot(gs[1])
ax2 = sns.barplot(hue="Condition", y="Signal.+", x="nodeClass", data=allCap, order=hueOrder, hue_order=["VEH", "DHT"], palette=colorPalette)
plt.xticks(rotation=45)
plt.legend()

fig.savefig(f"{figureRoot}/ProCap.bar.pdf")









































































filterCapG = capG[~((capG["SignalDHT"] == 0) & (capG["SignalVEH"] == 0))]

###########



fig = plt.figure(figsize=[15,6])
plt.subplots_adjust(bottom = 0.2)
gs = gridspec.GridSpec(ncols=2, nrows=1)


fig.add_subplot(gs[0])
sns.scatterplot(x="SignalVEH.+", y="SignalDHT.+", data=capG, hue="nodeClass", hue_order=hueOrder)
plt.legend()

fig.add_subplot(gs[1])
sns.scatterplot(x="SignalVEH.-", y="SignalDHT.-", data=capG, hue="nodeClass", hue_order=hueOrder)
plt.legend()


fig.savefig(f"{figureRoot}/ProCap.ARBS.Pro.scat.strand.pdf")


hueArbs = ["con", "ind", "non"]

fig = plt.figure(figsize=[15,6])
plt.subplots_adjust(bottom = 0.2)
gs = gridspec.GridSpec(ncols=2, nrows=1)


fig.add_subplot(gs[0])
sns.scatterplot(x="SignalVEH.+", y="SignalDHT.+", data=capG, hue="nodeClass", hue_order=hueArbs)
plt.legend()

fig.add_subplot(gs[1])
sns.scatterplot(x="SignalVEH.-", y="SignalDHT.-", data=capG, hue="nodeClass", hue_order=hueArbs)
plt.legend()


fig.savefig(f"{figureRoot}/ProCap.ARBS.scat.strand.pdf")

###### Without Strand

fig = plt.figure(figsize=[15,6])
plt.subplots_adjust(bottom = 0.2)
gs = gridspec.GridSpec(ncols=2, nrows=1)

fig.add_subplot(gs[0])
sns.scatterplot(x="SignalVEH", y="SignalDHT", data=capG, hue="nodeClass", hue_order=hueArbs)
plt.legend()

fig.add_subplot(gs[1])
sns.scatterplot(x="SignalVEH", y="SignalDHT", data=capG, hue="nodeClass", hue_order=hueOrder)
plt.legend()


fig.savefig(f"{figureRoot}/ProCap.ARBS.Pro.scat.pdf")















# fig = plt.figure(figsize=[9,12])
# ax = sns.boxplot(x="nodeClass", y="DHT", hue="Annotation", data=capG)
# ax.set_yscale("log")
# fig.savefig(f"{figureRoot}/ProCap.ARBS.pdf")


# fig = plt.figure(figsize=[15,9])
# ax = sns.barplot(hue="nodeClass", y="NumberofCapsDHT", x="Annotation", data=capG)
# # ax.set_yscale("log")
# fig.savefig(f"{figureRoot}/ProCap.hueARBS.bar.pdf")


hueOrder = ["con", "ind", "non", "tsP", "tsM"]
order = ["3' UTR", "5' UTR", "Distal Intergenic", "Downstream ",  "Exon ", "Intron ", "Promoter "]


fig = plt.figure(figsize=[15,6])
plt.subplots_adjust(bottom = 0.2)
gs = gridspec.GridSpec(ncols=2, nrows=1)
fig.add_subplot(gs[0])
ax1 = sns.barplot(hue="nodeClass", y="SignalVEH", x="annotation", data=filterCapG, order=order, hue_order=hueOrder)
plt.xticks(rotation=45)
plt.legend()

fig.add_subplot(gs[1])
ax2 = sns.barplot(hue="nodeClass", y="SignalDHT", x="annotation", data=filterCapG, order=order, hue_order=hueOrder)
plt.xticks(rotation=45)
plt.legend()

fig.savefig(f"{figureRoot}/ProCap.ARBS.Pro.bar.CS.filtered.pdf")



# boxplot


hueOrder = ["con", "ind", "non", "tsP", "tsM"]
order = ["3' UTR", "5' UTR", "Distal Intergenic", "Downstream ",  "Exon ", "Intron ", "Promoter "]
boxPairs = [(h, x) for x in order for h in hueOrder]

fig = plt.figure(figsize=[15,6])
plt.subplots_adjust(bottom = 0.2)
gs = gridspec.GridSpec(ncols=2, nrows=1)
fig.add_subplot(gs[0])
ax1 = sns.boxplot(hue="nodeClass", y="SignalVEH", x="annotation", data=filterCapG, order=order, hue_order=hueOrder)
plt.xticks(rotation=45)
plt.legend()
ax1.set_yscale("log")
# add_stat_annotation(ax1,x="annotation", y="SignalVEH", data=capG, hue="nodeClass",  test="Mann-Whitney", box_pairs=boxPairs, loc='inside', line_height=0)

fig.add_subplot(gs[1])
ax2 = sns.boxplot(hue="nodeClass", y="SignalDHT", x="annotation", data=filterCapG, order=order, hue_order=hueOrder)
plt.xticks(rotation=45)
plt.legend()
ax2.set_yscale("log")
# add_stat_annotation(ax2,x="annotation", y="SignalDHT", data=capG, hue="nodeClass",  test="Mann-Whitney", box_pairs=boxPairs, loc='inside', line_height=0)

fig.savefig(f"{figureRoot}/ProCap.ARBS.Pro.box.CS.pdf")








geneStrands = [1, 2]

titles = ["Plus Strand", "Minus Strand"]


fig = plt.figure(figsize=[30, 15])
plt.subplots_adjust(bottom = 0.2, wspace=0.3, hspace=0.4)
gs = gridspec.GridSpec(ncols=4, nrows=2)
for geneStrand in geneStrands:
    i = geneStrand -1
    fig.add_subplot(gs[0, i])
    ax1 = sns.barplot(hue="nodeClass", y="SignalVEH.+", x="annotation", data=capG[capG["geneStrand"] == geneStrand], hue_order=hueOrder, order=order)
    plt.xticks(rotation=45)
    plt.legend()
    plt.title(titles[i])
    plt.ylim([0,7000])
    fig.add_subplot(gs[1, i])
    ax2 = sns.barplot(hue="nodeClass", y="SignalDHT.+", x="annotation", data=capG[capG["geneStrand"] == geneStrand], hue_order=hueOrder, order=order)
    plt.xticks(rotation=45)
    plt.legend()
    plt.title(titles[i])
    plt.ylim([0,7000])
    fig.add_subplot(gs[0, i+2])
    ax1 = sns.barplot(hue="nodeClass", y="SignalVEH.-", x="annotation", data=capG[capG["geneStrand"] == geneStrand], hue_order=hueOrder, order=order)
    plt.xticks(rotation=45)
    plt.legend()
    plt.title(titles[i])
    plt.ylim([0,7000])
    fig.add_subplot(gs[1, i+2])
    ax2 = sns.barplot(hue="nodeClass", y="SignalDHT.-", x="annotation", data=capG[capG["geneStrand"] == geneStrand], hue_order=hueOrder, order=order)
    plt.xticks(rotation=45)
    plt.legend()
    plt.title(titles[i])
    plt.ylim([0,7000])

fig.savefig(f"{figureRoot}/ProCap.ARBS.Pro.bar.strand.CS.GS.pdf")



df3 = capG.groupby(["nodeClass", "annotation", "geneStrand"]).size().reset_index()
df4 = df3.pivot_table(index=["nodeClass", "geneStrand"], columns="annotation", values=0, aggfunc=np.sum)




capG = capG.sort_values("SYMBOL").reset_index(drop=True)

capG["qVal"] = 1
capG["LogFC"] = 0


capG["SYMBOL"] = capG["SYMBOL"].astype(str)
capG = capG.sort_values("SYMBOL").reset_index(drop=True)
logFC["Gene"] = logFC["Gene"].astype(str)
logFC = logFC .sort_values("Gene").reset_index(drop=True)

for i in range(len(capG)):
    sym = capG.loc[i, "SYMBOL"]
    leftL = bi.bisect_left(logFC["Gene"] , sym)
    rightL = bi.bisect_right(logFC["Gene"] , sym)
    tmpL = logFC[leftV:rightV]
    capG.loc[i, "qVal"] = logFC.loc[leftL:rightL, "padj"].mean()
    capG.loc[i, "LogFC"] = logFC.loc[leftL:rightL, "log2FoldChange"].mean()



capG["DE"] = "not"

capG.loc[(
(capG["LogFC"] > 1) &
(capG["qVal"] < 0.05)
), "DE"] = "upG"

capG.loc[(
(capG["LogFC"] < -1) &
(capG["qVal"] < 0.05)
), "DE"] = "dwG"




geneStrands = [1, 2]

titles = ["Plus Strand", "Minus Strand"]

DE = ["upG", "dwG"]

fig = plt.figure(figsize=[30, 30])
plt.subplots_adjust(bottom = 0.2, wspace=0.3, hspace=0.4)
gs = gridspec.GridSpec(ncols=4, nrows=4)

i = 0
for  geneStrand in geneStrands:
    for de in DE:
        fig.add_subplot(gs[i])
        ax1 = sns.barplot(hue="nodeClass", y="SignalVEH.+", x="annotation", data=capG[(capG["geneStrand"] == geneStrand) & (capG["DE"] == de) ], hue_order=["con", "ind", "non"], order=order)
        plt.xticks(rotation=45)
        plt.legend()
        plt.title(titles[geneStrand-1] + de)
        # plt.ylim([0,7000])
        i += 1
        fig.add_subplot(gs[i])
        ax2 = sns.barplot(hue="nodeClass", y="SignalDHT.+", x="annotation", data=capG[(capG["geneStrand"] == geneStrand) & (capG["DE"] == de)], hue_order=["con", "ind", "non"], order=order)
        plt.xticks(rotation=45)
        plt.legend()
        plt.title(titles[geneStrand-1] + de)
        # plt.ylim([0,7000])
        i += 1
        fig.add_subplot(gs[i])
        ax1 = sns.barplot(hue="nodeClass", y="SignalVEH.-", x="annotation", data=capG[(capG["geneStrand"] == geneStrand) & (capG["DE"] == de)], hue_order=["con", "ind", "non"], order=order)
        plt.xticks(rotation=45)
        plt.legend()
        plt.title(titles[geneStrand-1] + de)
        # plt.ylim([0,7000])
        i += 1
        fig.add_subplot(gs[i])
        ax2 = sns.barplot(hue="nodeClass", y="SignalDHT.-", x="annotation", data=capG[(capG["geneStrand"] == geneStrand) & (capG["DE"] == de)], hue_order=["con", "ind", "non"], order=order)
        plt.xticks(rotation=45)
        plt.legend()
        plt.title(titles[geneStrand-1] + de)
        # plt.ylim([0,7000])
        i += 1

fig.savefig(f"{figureRoot}/ProCap.ARBS.bar.strand.CS.GS.DE.pdf")




##################

capGdf = pd.DataFrame()

for i, edge in enumerate(G.edges()):
    aNode, bNode = edge
    achr, astart, aend = G.nodes[aNode]["nodeRange"]
    bchr, bstart, bend = G.nodes[bNode]["nodeRange"]
    capGdf.loc[i, "aChr"] = achr
    capGdf.loc[i, "aStart"] = astart
    capGdf.loc[i, "aEnd"] = aend
    capGdf.loc[i, "aName"] = aNode
    capGdf.loc[i, "anodeClass"] = G.nodes[aNode]["nodeClass"]
    capGdf.loc[i, "bChr"] = bchr
    capGdf.loc[i, "bStart"] = bstart
    capGdf.loc[i, "bEnd"] = bend
    capGdf.loc[i, "bName"] = bNode
    capGdf.loc[i, "bnodeClass"] = G.nodes[bNode]["nodeClass"]
    leftV = bi.bisect_left(veh["start.veh"] , astart)
    rightV = bi.bisect_left(veh["end.veh"] , aend)
    leftD = bi.bisect_left(dht["start.dht"] , astart)
    rightD = bi.bisect_left(dht["end.dht"] , aend)
    tmpV = veh[leftV:rightV]
    tmpV = tmpV[tmpV["chr.veh"] == achr]
    if tmpV.shape[0] > 0:
        capGdf.loc[i, "SignalVEH.+.a"] = max(max(list(tmpV[tmpV["strand.veh"] == "+"]["signal.veh"]), [0]))
        capGdf.loc[i, "SignalVEH.-.a"] = max(max(list(tmpV[tmpV["strand.veh"] == "-"]["signal.veh"]), [0]))
        capGdf.loc[i, "NumberofCapsVEH.+.a"] = tmpV[tmpV["strand.veh"] == "+"].shape[0]
        capGdf.loc[i, "NumberofCapsVEH.-.a"] = tmpV[tmpV["strand.veh"] == "-"].shape[0]
        capGdf.loc[i,"SignalVEH.a"] = tmpV["signal.veh"].max()
        capGdf.loc[i, "NumberofCapsVEH.a"] = tmpV.shape[0]
    tmpD = dht[leftD:rightD]
    tmpD = tmpD[tmpD["chr.dht"] == chr]
    if tmpD.shape[0] > 0:
        capGdf.loc[i, "SignalDHT.+.a"] = max(max(list(tmpD[tmpD["strand.dht"] == "+"]["signal.dht"]), [0]))
        capGdf.loc[i, "SignalDHT.-.a"] = max(max(list(tmpD[tmpD["strand.dht"] == "-"]["signal.dht"]), [0]))
        capGdf.loc[i, "NumberofCapsDHT.+.a"] = tmpD[tmpD["strand.dht"] == "+"].shape[0]
        capGdf.loc[i, "NumberofCapsDHT.-.a"] = tmpD[tmpD["strand.dht"] == "-"].shape[0]
        capGdf.loc[i,"SignalDHT.a"] = tmpD["signal.dht"].max()
        capGdf.loc[i, "NumberofCaplDHT.a"] = tmpD.shape[0]
    leftV = bi.bisect_left(veh["start.veh"] , bstart)
    rightV = bi.bisect_left(veh["end.veh"] , bend)
    leftD = bi.bisect_left(dht["start.dht"] , bstart)
    rightD = bi.bisect_left(dht["end.dht"] , bend)
    tmpV = veh[leftV:rightV]
    tmpV = tmpV[tmpV["chr.veh"] == bchr]
    if tmpV.shape[0] > 0:
        capGdf.loc[i, "SignalVEH.+.b"] = max(max(list(tmpV[tmpV["strand.veh"] == "+"]["signal.veh"]), [0]))
        capGdf.loc[i, "SignalVEH.-.b"] = max(max(list(tmpV[tmpV["strand.veh"] == "-"]["signal.veh"]), [0]))
        capGdf.loc[i, "NumberofCapsVEH.+.b"] = tmpV[tmpV["strand.veh"] == "+"].shape[0]
        capGdf.loc[i, "NumberofCapsVEH.-.b"] = tmpV[tmpV["strand.veh"] == "-"].shape[0]
        capGdf.loc[i,"SignalVEH.b"] = tmpV["signal.veh"].max()
        capGdf.loc[i, "NumberofCapsVEH.b"] = tmpV.shape[0]
    tmpD = dht[leftD:rightD]
    tmpD = tmpD[tmpD["chr.dht"] == chr]
    if tmpD.shape[0] > 0:
        capGdf.loc[i, "SignalDHT.+.b"] = max(max(list(tmpD[tmpD["strand.dht"] == "+"]["signal.dht"]), [0]))
        capGdf.loc[i, "SignalDHT.-.b"] = max(max(list(tmpD[tmpD["strand.dht"] == "-"]["signal.dht"]), [0]))
        capGdf.loc[i, "NumberofCapsDHT.+.b"] = tmpD[tmpD["strand.dht"] == "+"].shape[0]
        capGdf.loc[i, "NumberofCapsDHT.-.b"] = tmpD[tmpD["strand.dht"] == "-"].shape[0]
        capGdf.loc[i,"SignalDHT.b"] = tmpD["signal.dht"].max()
        capGdf.loc[i, "NumberofCaplDHT.b"] = tmpD.shape[0]


capGdf = capGdf.fillna(0)



aUp = capGdf[(capGdf["anodeClass"] == "upP")]
bUp = capGdf[(capGdf["bnodeClass"] == "upP")]

aUp["DHT"] = aUp["SignalDHT.b"] / aUp["SignalDHT.a"]+1
bUp["DHT"] = bUp["SignalDHT.a"] / bUp["SignalDHT.b"]+1


aUp["VEH"] = aUp["SignalVEH.b"] / aUp["SignalVEH.a"]+1
bUp["VEH"] = bUp["SignalVEH.a"] / bUp["SignalVEH.b"]+1


fig = plt.figure(figsize=[15, 15])
plt.subplots_adjust(bottom = 0.2, wspace=0.3, hspace=0.4)
gs = gridspec.GridSpec(ncols=2, nrows=2)
fig.add_subplot(gs[0])
sns.barplot(y="VEH", x="bnodeClass", data=aUp)
fig.add_subplot(gs[1])
sns.barplot(y="DHT", x="bnodeClass", data=aUp)
fig.add_subplot(gs[2])
sns.barplot(y="VEH", x="anodeClass", data=bUp)
fig.add_subplot(gs[3])
sns.barplot(y="DHT", x="anodeClass", data=bUp)
fig.savefig(f"{figureRoot}/ProCap.upInteraction.pdf")



# capG.groupby(["NumberofCapsDHT", "Annotation"]).size()
#
# fig = plt.figure(figsize=[15,6])
# plt.subplots_adjust(bottom = 0.2)
# gs = gridspec.GridSpec(ncols=2, nrows=1)
# fig.add_subplot(gs[0])
# ax1 = sns.barplot(hue="nodeClass", x="NumberofCapsVEH", x="Annotation", data=capG)
# plt.xticks(rotation=45)
# plt.legend()
#
# fig.add_subplot(gs[1])
# ax2 = sns.barplot(hue="nodeClass", x="NumberofCapsDHT", x="Annotation", data=capG)
# plt.xticks(rotation=45)
# plt.legend()
#
# fig.savefig(f"{figureRoot}/ProCap.ARBS.bar.Num.pdf")
