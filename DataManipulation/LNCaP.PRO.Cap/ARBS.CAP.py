
import pandas as pd
import bisect as bi
import numpy as np


home = "/home/ualtintas/LNCaP.PRO.cap"
arbsHome = "/home/ualtintas/ARBSs/regions"

# dataRoot = f"{home}/github/Data/CisGraph/Vers1.0"
# figureRoot = f"{home}/github/Figures/CisGraph/Vers1.0"

veh = pd.read_csv(f"{home}/GSM1543793_LNCaP_PRO-cap_1h_vehicle_hg18_peaks.bed", sep="\t", names=["chr", "start", "end", "score", "signal", "strand"])
dht = pd.read_csv(f"{home}/GSM1543794_LNCaP_PRO-cap_1h_dht_hg18_peaks.bed", sep="\t", names=["chr", "start", "end", "score", "signal", "strand"])

con = pd.read_csv(f"{arbsHome}/cons-arbs.bed", sep="\t", names=["chr", "start", "end", "name"])
ind = pd.read_csv(f"{arbsHome}/ind-arbs.bed", sep="\t", names=["chr", "start", "end", "name"])
non = pd.read_csv(f"{arbsHome}/Non-Active-ARBS.bed", sep="\t", names=["chr", "start", "end", "name"])
oth = pd.read_csv(f"{arbsHome}/otherARBS.bed", sep="\t", names=["chr", "start", "end", "name"])



veh = veh.sort_values("start").reset_index(drop=True)
dht = dht.sort_values("start").reset_index(drop=True)



con["caps.DHT"] = None
con["caps.VEH"] = None


for i in range(len(con)):
    cChr = con.loc[i, "chr"]
    cStart = con.loc[i, "start"] -n
    cEnd = con.loc[i, "end"] +n
    left = bi.bisect_left(veh["start"] , cStart)
    right = bi.bisect_left(veh["end"] , cEnd)
    tmp = veh[left:right]
    tmp = tmp[tmp["chr"]==cChr]
    con.loc[i,"caps.VEH"] = len(tmp)
    left = bi.bisect_left(dht["start"] , cStart)
    right = bi.bisect_left(dht["end"] , cEnd)
    tmp = dht[left:right]
    tmp = tmp[tmp["chr"]==cChr]
    con.loc[i, "caps.DHT"] = len(tmp)





ind["caps.DHT"] = None
ind["caps.VEH"] = None


for i in range(len(ind)):
    cChr = ind.loc[i, "chr"]
    cStart = ind.loc[i, "start"] -n
    cEnd = ind.loc[i, "end"] +n
    left = bi.bisect_left(veh["start"] , cStart)
    right = bi.bisect_left(veh["end"] , cEnd)
    tmp = veh[left:right]
    tmp = tmp[tmp["chr"]==cChr]
    ind.loc[i,"caps.VEH"] = len(tmp)
    left = bi.bisect_left(dht["start"] , cStart)
    right = bi.bisect_left(dht["end"] , cEnd)
    tmp = dht[left:right]
    tmp = tmp[tmp["chr"]==cChr]
    ind.loc[i, "caps.DHT"] = len(tmp)





non["caps.DHT"] = None
non["caps.VEH"] = None


for i in range(len(non)):
    cChr = non.loc[i, "chr"]
    cStart = non.loc[i, "start"] -n
    cEnd = non.loc[i, "end"] +n
    left = bi.bisect_left(veh["start"] , cStart)
    right = bi.bisect_left(veh["end"] , cEnd)
    tmp = veh[left:right]
    tmp = tmp[tmp["chr"]==cChr]
    non.loc[i,"caps.VEH"] = len(tmp)
    left = bi.bisect_left(dht["start"] , cStart)
    right = bi.bisect_left(dht["end"] , cEnd)
    tmp = dht[left:right]
    tmp = tmp[tmp["chr"]==cChr]
    non.loc[i, "caps.DHT"] = len(tmp)
