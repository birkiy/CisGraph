
from Functions.Packages import *
from Functions.PlotFunctions import *




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





veh = pd.read_csv(f"{dataRoot}/LNCaP.PRO.Cap/hglft_GSM1543793_LNCaP_PRO-cap_1h_vehicle_hg19_peaks.bed", sep="\t", names=["chr.veh", "start.veh", "end.veh", "score.veh", "signal.veh", "strand.veh"])
dht = pd.read_csv(f"{dataRoot}/LNCaP.PRO.Cap/hglft_GSM1543794_LNCaP_PRO-cap_1h_dht_hg19_peaks.bed", sep="\t", names=["chr.dht", "start.dht", "end.dht", "score.dht", "signal.dht", "strand.dht"])


veh = veh.sort_values("start.veh").reset_index(drop=True)
dht = dht.sort_values("start.dht").reset_index(drop=True)






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


tsPAnn = tsPAnn.reset_index(drop=True)

veh2 = pd.DataFrame()
dht2 = pd.DataFrame()

for i in range(len(tsPAnn)):
    chr = tsPAnn.loc[i, "Chr"]
    start = tsPAnn.loc[i, "Start"] + 250
    end = tsPAnn.loc[i, "End"] - 250
    leftV = bi.bisect_left(veh["start.veh"] , start)
    rightV = bi.bisect_left(veh["end.veh"] , end)
    leftD = bi.bisect_left(dht["start.dht"] , start)
    rightD = bi.bisect_left(dht["end.dht"] , end)
    tmpV = veh[leftV:rightV]
    tmpV = tmpV[tmpV["chr.veh"] == chr]
    veh2 = pd.concat([veh2, tmpV])
    if tmpV.shape[0] > 0:
        tsPAnn.loc[i, "SignalVEH.+"] = max(max(list(tmpV[tmpV["strand.veh"] == "+"]["signal.veh"]), [0]))
        tsPAnn.loc[i, "SignalVEH.-"] = max(max(list(tmpV[tmpV["strand.veh"] == "-"]["signal.veh"]), [0]))
        tsPAnn.loc[i, "NumberofCapsVEH.+"] = tmpV[tmpV["strand.veh"] == "+"].shape[0]
        tsPAnn.loc[i, "NumberofCapsVEH.-"] = tmpV[tmpV["strand.veh"] == "-"].shape[0]
        tsPAnn.loc[i,"SignalVEH"] = tmpV["signal.veh"].max()
        tsPAnn.loc[i, "NumberofCapsVEH"] = tmpV.shape[0]
    tmpD = dht[leftD:rightD]
    tmpD = tmpD[tmpD["chr.dht"] == chr]
    dht2 = pd.concat([dht2, tmpD])
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
    start = tsMAnn.loc[i, "Start"] + 250
    end = tsMAnn.loc[i, "End"] - 250
    leftV = bi.bisect_left(veh["start.veh"] , start)
    rightV = bi.bisect_left(veh["end.veh"] , end)
    leftD = bi.bisect_left(dht["start.dht"] , start)
    rightD = bi.bisect_left(dht["end.dht"] , end)
    tmpV = veh[leftV:rightV]
    tmpV = tmpV[tmpV["chr.veh"] == chr]
    veh2 = pd.concat([veh2, tmpV])
    if tmpV.shape[0] > 0:
        tsMAnn.loc[i, "SignalVEH.+"] = max(max(list(tmpV[tmpV["strand.veh"] == "+"]["signal.veh"]), [0]))
        tsMAnn.loc[i, "SignalVEH.-"] = max(max(list(tmpV[tmpV["strand.veh"] == "-"]["signal.veh"]), [0]))
        tsMAnn.loc[i, "NumberofCapsVEH.+"] = tmpV[tmpV["strand.veh"] == "+"].shape[0]
        tsMAnn.loc[i, "NumberofCapsVEH.-"] = tmpV[tmpV["strand.veh"] == "-"].shape[0]
        tsMAnn.loc[i,"SignalVEH"] = tmpV["signal.veh"].max()
        tsMAnn.loc[i, "NumberofCapsVEH"] = tmpV.shape[0]
    tmpD = dht[leftD:rightD]
    tmpD = tmpD[tmpD["chr.dht"] == chr]
    dht2 = pd.concat([dht2, tmpD])
    if tmpD.shape[0] > 0:
        tsMAnn.loc[i, "SignalDHT.+"] = max(max(list(tmpD[tmpD["strand.dht"] == "+"]["signal.dht"]), [0]))
        tsMAnn.loc[i, "SignalDHT.-"] = max(max(list(tmpD[tmpD["strand.dht"] == "-"]["signal.dht"]), [0]))
        tsMAnn.loc[i, "NumberofCapsDHT.+"] = tmpD[tmpD["strand.dht"] == "+"].shape[0]
        tsMAnn.loc[i, "NumberofCapsDHT.-"] = tmpD[tmpD["strand.dht"] == "-"].shape[0]
        tsMAnn.loc[i,"SignalDHT"] = tmpD["signal.dht"].max()
        tsMAnn.loc[i, "NumberofCaplDHT"] = tmpD.shape[0]





veh2 = veh2.drop_duplicates()
dht2 = dht2.drop_duplicates()

veh3 = veh.drop(veh2.index, axis=0)
dht3 = dht.drop(dht2.index, axis=0)




veh3 = veh.reset_index(drop=True)
dht3 = dht.reset_index(drop=True)









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









con = conAnn[["V4", "Chr", "Start", "End", "annotation", "geneStrand", "SignalVEH.-", "NumberofCapsVEH.-", "SignalVEH.+", "NumberofCapsVEH.+", "SignalVEH", "NumberofCapsVEH", "SignalDHT.-", "NumberofCapsDHT.-", "SignalDHT.+", "NumberofCapsDHT.+", "SignalDHT", "NumberofCapsDHT", "SYMBOL"]]
ind = indAnn[["V4", "Chr", "Start", "End", "annotation", "geneStrand", "SignalVEH.-", "NumberofCapsVEH.-", "SignalVEH.+", "NumberofCapsVEH.+", "SignalVEH", "NumberofCapsVEH", "SignalDHT.-", "NumberofCapsDHT.-", "SignalDHT.+", "NumberofCapsDHT.+", "SignalDHT", "NumberofCapsDHT", "SYMBOL"]]
non = nonAnn[["V4", "Chr", "Start", "End", "annotation", "geneStrand", "SignalVEH.-", "NumberofCapsVEH.-", "SignalVEH.+", "NumberofCapsVEH.+", "SignalVEH", "NumberofCapsVEH", "SignalDHT.-", "NumberofCapsDHT.-", "SignalDHT.+", "NumberofCapsDHT.+", "SignalDHT", "NumberofCapsDHT", "SYMBOL"]]
# tsP = tsPAnn[["V4", "Chr", "Start", "End", "annotation", "geneStrand", "SignalVEH.-", "NumberofCapsVEH.-", "SignalVEH.+", "NumberofCapsVEH.+", "SignalVEH", "NumberofCapsVEH", "SignalDHT.-", "NumberofCapsDHT.-", "SignalDHT.+", "NumberofCapsDHT.+", "SignalDHT", "NumberofCapsDHT", "SYMBOL"]]
# tsM = tsMAnn[["V4", "Chr", "Start", "End", "annotation", "geneStrand", "SignalVEH.-", "NumberofCapsVEH.-", "SignalVEH.+", "NumberofCapsVEH.+", "SignalVEH", "NumberofCapsVEH", "SignalDHT.-", "NumberofCapsDHT.-", "SignalDHT.+", "NumberofCapsDHT.+", "SignalDHT", "NumberofCapsDHT", "SYMBOL"]]


con["nodeClass"] = "con"
ind["nodeClass"] = "ind"
non["nodeClass"] = "non"
# tsP["nodeClass"] = "tsP"
# tsM["nodeClass"] = "tsM"


# capG = pd.concat([con, ind, non, tsP, tsM]).reset_index(drop=True)
capG = pd.concat([con, ind, non]).reset_index(drop=True)

capG["annotation"] = capG["annotation"].str.split("(", expand=True)[0]

df = capG.groupby(["nodeClass", "annotation"]).size().reset_index()

df.pivot(index="annotation", columns="nodeClass", values=0)








capSig = capG[["nodeClass", "annotation", "SignalDHT.+", "SignalDHT.-", "SignalVEH.+", "SignalVEH.-", "V4"]]


dhtCap = capSig[["nodeClass", "annotation", "SignalDHT.+", "SignalDHT.-", "V4"]]
dhtCap.loc[:, "Condition"] = "DHT"
dhtCap = dhtCap.rename(columns={"SignalDHT.+": "Signal.+", "SignalDHT.-": "Signal.-"})

vehCap = capSig[["nodeClass", "annotation", "SignalVEH.+", "SignalVEH.-", "V4"]]
vehCap.loc[:,"Condition"] = "VEH"
vehCap = vehCap.rename(columns={"SignalVEH.+": "Signal.+", "SignalVEH.-": "Signal.-"})

allCap = pd.concat([dhtCap,vehCap])


allCap["mean"] = allCap["Signal.+"] + allCap["Signal.-"]


nodeClass = ["con", "ind", "non"]

filtCap = allCap[allCap["mean"] != 0]


for node1 in nodeClass:
    print(node1)
    scipy.stats.ranksums(filtCap.loc[(filtCap["nodeClass"] == node1) & (filtCap["Condition"] == "VEH"), "mean"], filtCap.loc[(filtCap["nodeClass"] == node1) & (filtCap["Condition"] == "DHT"), "mean"])
    for node2 in nodeClass:
        print(node1, node2)
        scipy.stats.ranksums(filtCap.loc[filtCap["nodeClass"] == node1, "mean"], filtCap.loc[filtCap["nodeClass"] == node2, "mean"])
        scipy.stats.ranksums(filtCap.loc[filtCap["nodeClass"] == node1, "mean"], filtCap.loc[filtCap["nodeClass"] == node2, "mean"])






allCap["Signal.+"] += 1
allCap["Signal.-"] += 1


allCap["Signal.+"] -= 1
allCap["Signal.-"] -= 1



allCap["annotation2"] = "Others"
allCap.loc[allCap["annotation"] == "Promoter ", "annotation2"] = "Promoter"
allCap["annotation3"] = allCap[["annotation2", "Condition"]].agg(".".join, axis=1)



fig = plt.figure(figsize=[15,6])
plt.subplots_adjust(bottom = 0.2)
gs = gridspec.GridSpec(ncols=2, nrows=1)


params = dict(data=allCap,
              x='nodeClass',
              hue='Condition',
              dodge=True,
              order=["con", "ind", "non"],
              hue_order=["VEH", "DHT"])

ax1 = fig.add_subplot(gs[0])
ax1 = sns.stripplot(y="Signal.-", **params)
# ax1Box = sns.boxplot(y="Signal.-", palette=['#BBBBBB','#DDDDDD'], **params)
ax1.set(yscale="log")
ax1.set_ylim([1,2*10**5])
plt.legend()


ax2 = fig.add_subplot(gs[1])
ax2 = sns.stripplot(y="Signal.+", **params)
# ax2Box = sns.boxplot(y="Signal.+", palette=['#BBBBBB','#DDDDDD'], **params)
ax2.set(yscale="log")
ax2.set_ylim([1,2*10**5])
plt.legend()

fig.savefig(f"{figureRoot}/Pro.NewStrip.pdf")


###########################################################################



fig = plt.figure(figsize=[15,6])
plt.subplots_adjust(bottom = 0.2)
gs = gridspec.GridSpec(ncols=2, nrows=1)


params = dict(data=allCap,
              x='nodeClass',
              hue='annotation3',
              dodge=True,
              order=["con", "ind", "non"],
              palette="Set2")

ax1 = fig.add_subplot(gs[0])
ax1 = sns.stripplot(y="Signal.-", **params)
# ax1Box = sns.boxplot(y="Signal.-", palette=['#BBBBBB','#DDDDDD'], **params)
ax1.set(yscale="log")
ax1.set_ylim([1,2*10**5])



ax2 = fig.add_subplot(gs[1])
ax2 = sns.stripplot(y="Signal.+", **params)
# ax2Box = sns.boxplot(y="Signal.+", palette=['#BBBBBB','#DDDDDD'], **params)
ax2.set(yscale="log")
ax2.set_ylim([1,2*10**5])
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax1.get_legend().remove()

fig.savefig(f"{figureRoot}/Pro.NewStrip.onto.pdf")


###########################################################################


nodeClass = ["con", "ind", "non", "tsP", "tsM"]

for node in nodeClass:
    dhtP = dhtCap[dhtCap["nodeClass"] == node]["Signal.+"]
    dhtM = dhtCap[dhtCap["nodeClass"] == node]["Signal.-"]
    vehP = vehCap[vehCap["nodeClass"] == node]["Signal.+"]
    vehM = vehCap[vehCap["nodeClass"] == node]["Signal.-"]
    s, p = scipy.stats.ranksums(dhtP, vehP)
    print(f'{node} class DHT vs VEH for + strand, statistic:{s}, pvalue:{p}')
    s, p = scipy.stats.ranksums(dhtM, vehM)
    print(f'{node} class DHT vs VEH for - strand, statistic:{s}, pvalue:{p}')



##############################################


allCap["Signal.+"] += 1
allCap["Signal.-"] += 1


fig = plt.figure(figsize=[15,6])
plt.subplots_adjust(bottom = 0.2)
gs = gridspec.GridSpec(ncols=2, nrows=1)


params = dict(data=allCap,
              x='nodeClass',
              hue='annotation3',
              dodge=True,
              order=["con", "ind", "non"],
              palette="Set2")

ax1 = fig.add_subplot(gs[0])
ax1 = sns.violinplot(y="Signal.-", **params)
# ax1Box = sns.boxplot(y="Signal.-", palette=['#BBBBBB','#DDDDDD'], **params)
ax1.set(yscale="log")
ax1.set_ylim([1,2*10**5])



ax2 = fig.add_subplot(gs[1])
ax2 = sns.violinplot(y="Signal.+", **params)
# ax2Box = sns.boxplot(y="Signal.+", palette=['#BBBBBB','#DDDDDD'], **params)
ax2.set(yscale="log")
ax2.set_ylim([1,2*10**5])
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax1.get_legend().remove()

fig.savefig(f"{figureRoot}/Pro.violin.pdf")







conGroDht = allDF[(allDF["class"] == "con") & (allDF["Condition"] == "DHT")]
conProDht = allCap[(allCap["nodeClass"] == "con") & (allCap["Condition"] == "DHT")]


conGroDht = conGroDht.sort_values("Gene")
conProDht = conProDht.sort_values("V4")



scipy.stats.spearmanr(conGroDht["GroSeqMinus"], conProDht["Signal.+"])
scipy.stats.spearmanr(conGroDht["GroSeqPlus"], conProDht["Signal.-"])







































allCap2 = allCap
allCap = vehCap

allCap = dhtCap

sum(allCap.loc[(allCap["nodeClass"] == "con"), "Signal.+"] > 0)
sum(allCap.loc[((allCap["nodeClass"] == "con") & (allCap["annotation"] == "Promoter ")), "Signal.+"] > 0)


sum(allCap.loc[(allCap["nodeClass"] == "ind"), "Signal.+"] > 0)
sum(allCap.loc[((allCap["nodeClass"] == "ind") & (allCap["annotation"] == "Promoter ")), "Signal.+"] > 0)


sum(allCap.loc[(allCap["nodeClass"] == "non"), "Signal.+"] > 0)
sum(allCap.loc[((allCap["nodeClass"] == "non") & (allCap["annotation"] == "Promoter ")), "Signal.+"] > 0)

# Signal -

sum(allCap.loc[(allCap["nodeClass"] == "con"), "Signal.-"] > 0)
sum(allCap.loc[((allCap["nodeClass"] == "con") & (allCap["annotation"] == "Promoter ")), "Signal.-"] > 0)


sum(allCap.loc[(allCap["nodeClass"] == "ind"), "Signal.-"] > 0)
sum(allCap.loc[((allCap["nodeClass"] == "ind") & (allCap["annotation"] == "Promoter ")), "Signal.-"] > 0)


sum(allCap.loc[(allCap["nodeClass"] == "non"), "Signal.-"] > 0)
sum(allCap.loc[((allCap["nodeClass"] == "non") & (allCap["annotation"] == "Promoter ")), "Signal.-"] > 0)



# Divergents

sum((allCap.loc[(allCap["nodeClass"] == "con"), "Signal.+"] > 0) & (allCap.loc[(allCap["nodeClass"] == "con"), "Signal.-"] > 0))

sum((allCap.loc[(allCap["nodeClass"] == "ind"), "Signal.+"] > 0) & (allCap.loc[(allCap["nodeClass"] == "ind"), "Signal.-"] > 0))

sum((allCap.loc[(allCap["nodeClass"] == "non"), "Signal.+"] > 0) & (allCap.loc[(allCap["nodeClass"] == "non"), "Signal.-"] > 0))



sum((allCap.loc[((allCap["nodeClass"] == "con") & (allCap["annotation"] == "Promoter ")), "Signal.+"] > 0) & (allCap.loc[((allCap["nodeClass"] == "con") & (allCap["annotation"] == "Promoter ")), "Signal.-"] > 0))

sum((allCap.loc[((allCap["nodeClass"] == "ind") & (allCap["annotation"] == "Promoter ")), "Signal.+"] > 0) & (allCap.loc[((allCap["nodeClass"] == "non") & (allCap["annotation"] == "Promoter ")), "Signal.-"] > 0))

sum((allCap.loc[((allCap["nodeClass"] == "non") & (allCap["annotation"] == "Promoter ")), "Signal.+"] > 0) & (allCap.loc[((allCap["nodeClass"] == "non") & (allCap["annotation"] == "Promoter ")), "Signal.-"] > 0))



# Describe


allCap.loc[(allCap["nodeClass"] == "con"), ["Signal.+", "Signal.-"]].describe()
allCap.loc[(allCap["nodeClass"] == "ind"), ["Signal.+", "Signal.-"]].describe()
allCap.loc[(allCap["nodeClass"] == "non"), ["Signal.+", "Signal.-"]].describe()

allCap.loc[((allCap["nodeClass"] == "con") & (allCap["annotation"] == "Promoter ")), ["Signal.+", "Signal.-"]].describe()
allCap.loc[((allCap["nodeClass"] == "ind") & (allCap["annotation"] == "Promoter ")), ["Signal.+", "Signal.-"]].describe()
allCap.loc[((allCap["nodeClass"] == "non") & (allCap["annotation"] == "Promoter ")), ["Signal.+", "Signal.-"]].describe()
