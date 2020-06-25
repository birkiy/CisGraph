

from Functions.Packages import *
from Functions.PlotFunctions import *

G = pickle.load(open(f"{dataRoot}/tmpData/GraphsGData.p", "rb" ))


veh = pd.read_csv(f"{dataRoot}/LNCaP.PRO.cap/GSM1543793_LNCaP_PRO-cap_1h_vehicle_hg18_peaks.bed", sep="\t", names=["chr.veh", "start.veh", "end.veh", "score.veh", "signal.veh", "strand.veh"])
dht = pd.read_csv(f"{dataRoot}/LNCaP.PRO.cap/GSM1543794_LNCaP_PRO-cap_1h_dht_hg18_peaks.bed", sep="\t", names=["chr.dht", "start.dht", "end.dht", "score.dht", "signal.dht", "strand.dht"])

conAnn = pd.read_csv("/home/ualtintas/ARBSs/annotations/con.Annotations.txt", sep="\t")
conAnn = conAnn.rename(columns={"PeakID (cmd=annotatePeaks.pl regions/cons-arbs.bed hg19)": "PeakID"})
indAnn = pd.read_csv("/home/ualtintas/ARBSs/annotations/ind.Annotations.txt", sep="\t")
indAnn = indAnn.rename(columns={"PeakID (cmd=annotatePeaks.pl regions/ind-arbs.bed hg19)": "PeakID"})
nonAnn = pd.read_csv("/home/ualtintas/ARBSs/annotations/non.Annotations.txt", sep="\t")
nonAnn = nonAnn.rename(columns={"PeakID (cmd=annotatePeaks.pl regions/Non-Active-ARBS.bed hg19)": "PeakID"})


veh = veh.sort_values("start.veh").reset_index(drop=True)
dht = dht.sort_values("start.dht").reset_index(drop=True)


capG = pd.DataFrame()
for node in G.nodes():
    chr, start, end = G.nodes[node]["nodeRange"]
    nc = G.nodes[node]["nodeClass"]
    leftV = bi.bisect_left(veh["start"] , start)
    rightV = bi.bisect_left(veh["end"] , end)
    leftD = bi.bisect_left(dht["start"] , start)
    rightD = bi.bisect_left(dht["end"] , end)
    tmpV = veh[leftV:rightV]
    tmpV = tmpV[tmpV["chr"] == chr]
    tmpD = dht[leftD:rightD]
    tmpD = tmpD[tmpD["chr"] == chr]
    pd.concat[]


conAnn["VEH"] = None
conAnn["DHT"] = None

capG = pd.DataFrame()
for i in range(len(conAnn)):
    chr = conAnn.loc[i, "Chr"]
    start = conAnn.loc[i, "Start"]
    end = conAnn.loc[i, "End"]
    leftV = bi.bisect_left(veh["start.veh"] , start)
    rightV = bi.bisect_left(veh["end.veh"] , end)
    leftD = bi.bisect_left(dht["start.dht"] , start)
    rightD = bi.bisect_left(dht["end.dht"] , end)
    tmpV = veh[leftV:rightV]
    tmpV = tmpV[tmpV["chr.veh"] == chr]
    if tmpV.shape[0] > 0:
        conAnn.loc[i, "VEH"] = max(tmpV["signal.veh"])
    tmpD = dht[leftD:rightD]
    tmpD = tmpD[tmpD["chr.dht"] == chr]
    if tmpD.shape[0] > 0:
        conAnn.loc[i, "DHT"] = max(tmpD["signal.dht"])


indAnn["VEH"] = None
indAnn["DHT"] = None


for i in range(len(indAnn)):
    chr = indAnn.loc[i, "Chr"]
    start = indAnn.loc[i, "Start"]
    end = indAnn.loc[i, "End"]
    leftV = bi.bisect_left(veh["start.veh"] , start)
    rightV = bi.bisect_left(veh["end.veh"] , end)
    leftD = bi.bisect_left(dht["start.dht"] , start)
    rightD = bi.bisect_left(dht["end.dht"] , end)
    tmpV = veh[leftV:rightV]
    tmpV = tmpV[tmpV["chr.veh"] == chr]
    if tmpV.shape[0] > 0:
        indAnn.loc[i, "VEH"] = max(tmpV["signal.veh"])
    tmpD = dht[leftD:rightD]
    tmpD = tmpD[tmpD["chr.dht"] == chr]
    if tmpD.shape[0] > 0:
        indAnn.loc[i, "DHT"] = max(tmpD["signal.dht"])



nonAnn["VEH"] = None
nonAnn["DHT"] = None


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
        nonAnn.loc[i, "VEH"] = max(tmpV["signal.veh"])
    tmpD = dht[leftD:rightD]
    tmpD = tmpD[tmpD["chr.dht"] == chr]
    if tmpD.shape[0] > 0:
        nonAnn.loc[i, "DHT"] = max(tmpD["signal.dht"])


con = conAnn[["PeakID", "Chr", "Start", "End", "Annotation", "Gene Type", "VEH", "DHT"]]
ind = indAnn[["PeakID", "Chr", "Start", "End", "Annotation", "Gene Type", "VEH", "DHT"]]
non = nonAnn[["PeakID", "Chr", "Start", "End", "Annotation", "Gene Type", "VEH", "DHT"]]


con["nodeClass"] = "con"
ind["nodeClass"] = "ind"
non["nodeClass"] = "non"


capG = pd.concat([con, ind, non]).reset_index(drop=True)

capG["Annotation"] = capG["Annotation"].str.split(" ", expand=True)[0]

fig = plt.figure(figsize=[9,12])
sns.boxplot(x="nodeClass", y="DHT", hue="Annotation", data=capG)

fig.savefig(f"{figureRoot}/ProCap.ARBS.pdf")
