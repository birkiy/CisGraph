


from Functions.Packages import *

chipARBS = pd.read_csv(f"{dataRoot}/ARBS_chip_seq_2.csv")

chipProm = pd.read_csv(f"{dataRoot}/Promoters_chip_seq.csv")

chip = pd.concat([chipARBS, chipProm])




conBed = pd.read_csv(home + "/ARBSs/regions/cons-arbs.bed", sep="\t", names=["chr", "start", "end", "name"])
indBed = pd.read_csv(home + "/ARBSs/regions/ind-arbs.bed", sep="\t", names=["chr", "start", "end", "name"])
nonBed = pd.read_csv(home + "/ARBSs/regions/Non-Active-ARBS.bed", sep="\t", names=["chr", "start", "end", "name"])
nARBed = pd.read_csv(home + "/ARBSs/regions/negativeControl.ARBS.bed", sep="\t", names=["chr", "start", "end", "name"])

tsPBed = pd.read_csv(home + "/genomeAnnotations/Regions/TSS.hg19.+.bed", sep="\t", names=["chr", "start", "end", "name", "strand"])
tsMBed = pd.read_csv(home + "/genomeAnnotations/Regions/TSS.hg19.-.bed", sep="\t", names=["chr", "start", "end", "name", "strand"])

conChip = chip[chip["Unnamed: 0"].isin(conBed["name"])]
conChip["nodeClass"] = "con"
indChip = chip[chip["Unnamed: 0"].isin(indBed["name"])]
indChip["nodeClass"] = "ind"
nonChip = chip[chip["Unnamed: 0"].isin(nonBed["name"])]
nonChip["nodeClass"] = "non"
nARChip = chip[chip["Unnamed: 0"].isin(nARBed["name"])]
nARChip["nodeClass"] = "nAR"

tsPChip = chip[chip["Unnamed: 0"].isin(tsPBed["name"])]
tsPChip["nodeClass"] = "pro"
tsMChip = chip[chip["Unnamed: 0"].isin(tsMBed["name"])]
tsMChip["nodeClass"] = "pro"

chipAll = pd.concat([conChip, indChip, nonChip, nARChip, tsPChip, tsMChip])

chipAll = chipAll.reset_index(drop=True)

pickle.dump(chipAll, open(f"{dataRoot}/ChipAll.p", "wb" ))

###########################

rnaPol = chipAll[["Unnamed: 0", "nodeClass", 'POLR2A-ETOH-GSE43253', "POLR2A-DHT-GSE43253"]]


etMe = rnaPol[["Unnamed: 0", "nodeClass", 'POLR2A-ETOH-GSE43253']]
etMe["Condition"] = "EtOH"
etMe = etMe.rename(columns={'POLR2A-ETOH-GSE43253': 'POLR2A-GSE43253'})



r1Me = rnaPol[["Unnamed: 0", "nodeClass", "POLR2A-DHT-GSE43253"]]
r1Me["Condition"] = "DHT"
r1Me = r1Me.rename(columns={'POLR2A-DHT-GSE43253': 'POLR2A-GSE43253'})


dfMe = pd.concat([etMe, r1Me])

#colorPalette = ["#63b7af", "#abf0e9", "#d4f3ef", "#f5fffd", "#ee8572"]
colorPalette = {"con": "#5A5A5A",
                "ind": "#F9746D",
                "non": "#ACACAC",
                "nAR": "#F5F5F5",
                "pro": "#000000"}


fig = plt.figure(figsize=[8, 5])
gs = gridspec.GridSpec(ncols=2, nrows=1)
plt.subplots_adjust(wspace=0.4)


params = dict(x='nodeClass',
              y='POLR2A-GSE43253',
              order=["con", "ind", "non", "nAR", "pro"])

boxPairs = [("nAR", "con"),
            ("con", "ind"),
            ("con", "non"),
            ("ind", "non"),
            ("con", "pro")
            ]

ax1 = fig.add_subplot(gs[0])
ax1 = sns.boxplot(**params, data=dfMe[dfMe["Condition"] == "EtOH"],  palette=colorPalette)
plt.title("EtOH")
plt.ylabel("RNAPII - Occupancy Score", fontsize=14)
plt.xlabel("Node Classes", fontsize=14)
for tick in ax1.xaxis.get_major_ticks():
    tick.label.set_fontsize(13)


add_stat_annotation(ax1, **params, data=dfMe[dfMe["Condition"] == "EtOH"], box_pairs=boxPairs, test='Mann-Whitney')

ax2 = fig.add_subplot(gs[1])
ax2 = sns.boxplot(**params, data=dfMe[dfMe["Condition"] == "DHT"], palette=colorPalette)
plt.title("DHT")
plt.ylabel("RNAPII - Occupancy Score", fontsize=14)
plt.xlabel("Node Classes", fontsize=14)
for tick in ax2.xaxis.get_major_ticks():
    tick.label.set_fontsize(13)


add_stat_annotation(ax2, **params, data=dfMe[dfMe["Condition"] == "DHT"], box_pairs=boxPairs, test='Mann-Whitney')


fig.savefig(f"{figureRoot}/RNAPII.pdf")


##########################################3


hMeMark = chipAll[["Unnamed: 0", "nodeClass", "H3K4me3-ETOH-GSE114732", "H3K4me3-R1881-GSE114732", "H3K4me1-ETOH-GSE114732", "H3K4me1-R1881-GSE114732"]]

hMeMark["H3K4me3/H3K4me1-ETOH"] = (hMeMark["H3K4me3-ETOH-GSE114732"]) / (hMeMark["H3K4me1-ETOH-GSE114732"] )
hMeMark["H3K4me3/H3K4me1-R1881"] = (hMeMark["H3K4me3-R1881-GSE114732"]) / (hMeMark["H3K4me1-R1881-GSE114732"])

etMe = hMeMark[["Unnamed: 0", "nodeClass", "H3K4me3/H3K4me1-ETOH", "H3K4me3-ETOH-GSE114732", "H3K4me1-ETOH-GSE114732"]]
etMe["Condition"] = "EtOH"
etMe = etMe.rename(columns={"H3K4me3/H3K4me1-ETOH": "H3K4me3/H3K4me1", "H3K4me3-ETOH-GSE114732": "H3K4me3-GSE114732", "H3K4me1-ETOH-GSE114732": "H3K4me1-GSE114732"})



r1Me = hMeMark[["Unnamed: 0", "nodeClass", "H3K4me3/H3K4me1-R1881", "H3K4me3-R1881-GSE114732", "H3K4me1-R1881-GSE114732"]]
r1Me["Condition"] = "R1881"
r1Me = r1Me.rename(columns={"H3K4me3/H3K4me1-R1881": "H3K4me3/H3K4me1", "H3K4me3-R1881-GSE114732": "H3K4me3-GSE114732", "H3K4me1-R1881-GSE114732": "H3K4me1-GSE114732"})


dfMe = pd.concat([etMe, r1Me])


dfMe['H3K4me3/H3K4me1'] = dfMe['H3K4me3/H3K4me1'].apply(lambda x: 0 if x == float("inf") else x)



fig = plt.figure(figsize=[8, 5])
gs = gridspec.GridSpec(ncols=2, nrows=1)
plt.subplots_adjust(wspace=0.4, hspace=0.3)

#colorPalette = ["#63b7af", "#abf0e9", "#d4f3ef", "#f5fffd", "#ee8572"]


boxPairs = [("nAR", "con"),
            ("con", "ind"),
            ("con", "non"),
            ("ind", "non"),
            ("con", "pro")
            ]

params = dict(x='nodeClass',
              order=["con", "ind", "non", "nAR", "pro"])


ax1 = fig.add_subplot(gs[0])
ax1 = sns.boxplot(data=dfMe[dfMe["Condition"] == "EtOH"], y="H3K4me3-GSE114732", palette=colorPalette, **params)
plt.title("EtOH")
plt.ylabel("H3K4me3 - Occupancy Score", fontsize=14)
plt.xlabel("Node Classes", fontsize=14)
for tick in ax1.xaxis.get_major_ticks():
    tick.label.set_fontsize(13)


add_stat_annotation(ax1, data=dfMe[dfMe["Condition"] == "EtOH"], y="H3K4me3-GSE114732", **params, box_pairs=boxPairs, test='Mann-Whitney')

ax2 = fig.add_subplot(gs[1])
ax2 = sns.boxplot(data=dfMe[dfMe["Condition"] == "R1881"], y="H3K4me3-GSE114732", palette=colorPalette, **params)
plt.title("R1881")
plt.ylabel("H3K4me3 - Occupancy Score", fontsize=14)
plt.xlabel("Node Classes", fontsize=14)
for tick in ax2.xaxis.get_major_ticks():
    tick.label.set_fontsize(13)


add_stat_annotation(ax2,data=dfMe[dfMe["Condition"] == "R1881"], y="H3K4me3-GSE114732", **params, box_pairs=boxPairs, test='Mann-Whitney')


fig.savefig(f"{figureRoot}/H3K4me3.ARBS.pdf")

fig = plt.figure(figsize=[8, 5])
gs = gridspec.GridSpec(ncols=2, nrows=1)
plt.subplots_adjust(wspace=0.4, hspace=0.3)


boxPairs += [("nAR", "pro")]

ax1 = fig.add_subplot(gs[0])
ax1 = sns.boxplot(data=dfMe[dfMe["Condition"] == "EtOH"], y="H3K4me1-GSE114732", palette=colorPalette, **params)
plt.title("EtOH")
plt.ylabel("H3K4me1 - Occupancy Score", fontsize=14)
plt.xlabel("Node Classes", fontsize=14)
for tick in ax1.xaxis.get_major_ticks():
    tick.label.set_fontsize(13)


add_stat_annotation(ax1, data=dfMe[dfMe["Condition"] == "EtOH"], y="H3K4me1-GSE114732", **params, box_pairs=boxPairs, test='Mann-Whitney')

ax2 = fig.add_subplot(gs[1])
ax2 = sns.boxplot(data=dfMe[dfMe["Condition"] == "R1881"], y="H3K4me1-GSE114732", palette=colorPalette, **params)
plt.title("R1881")
plt.ylabel("H3K4me1 - Occupancy Score", fontsize=14)
plt.xlabel("Node Classes", fontsize=14)
for tick in ax2.xaxis.get_major_ticks():
    tick.label.set_fontsize(13)

add_stat_annotation(ax2,data=dfMe[dfMe["Condition"] == "R1881"], y="H3K4me1-GSE114732", **params, box_pairs=boxPairs, test='Mann-Whitney')


fig.savefig(f"{figureRoot}/H3K4me1.ARBS.pdf")

##########################




hMeMark = chipAll[["Unnamed: 0", "nodeClass", "H3K27ac-DHT-GSE51621", "H3K27ac-ETOH-GSE51621"]]

etMe = hMeMark[["Unnamed: 0", "nodeClass", "H3K27ac-ETOH-GSE51621"]]
etMe["Condition"] = "EtOH"
etMe = etMe.rename(columns={"H3K27ac-ETOH-GSE51621": "H3K27ac"})



r1Me = hMeMark[["Unnamed: 0", "nodeClass", "H3K27ac-DHT-GSE51621"]]
r1Me["Condition"] = "DHT"
r1Me = r1Me.rename(columns={"H3K27ac-DHT-GSE51621": "H3K27ac"})


dfMe = pd.concat([etMe, r1Me])



fig = plt.figure(figsize=[8, 5])
gs = gridspec.GridSpec(ncols=2, nrows=1)
plt.subplots_adjust(wspace=0.4)

#colorPalette = ["#63b7af", "#abf0e9", "#d4f3ef", "#f5fffd", "#ee8572"]

params = dict(x='nodeClass',
              y='H3K27ac',
              order=["con", "ind", "non", "nAR", "pro"])

boxPairs = [("nAR", "con"),
            ("con", "ind"),
            ("con", "non"),
            ("ind", "non"),
            ("con", "pro"),
            ("non", "pro"),
            ("nAR", "pro")
            ]

ax1 = fig.add_subplot(gs[0])
ax1 = sns.boxplot(**params, data=dfMe[dfMe["Condition"] == "EtOH"], palette=colorPalette)
plt.title("EtOH")
plt.ylabel("H3K27ac - Occupancy Score", fontsize=14)
plt.xlabel("Node Classes", fontsize=14)
for tick in ax1.xaxis.get_major_ticks():
    tick.label.set_fontsize(13)


add_stat_annotation(ax1, **params, data=dfMe[dfMe["Condition"] == "EtOH"], box_pairs=boxPairs, test='Mann-Whitney')

ax2 = fig.add_subplot(gs[1])
ax2 = sns.boxplot(**params, data=dfMe[dfMe["Condition"] == "DHT"], palette=colorPalette)
plt.title("DHT")
plt.ylabel("H3K27ac - Occupancy Score", fontsize=14)
plt.xlabel("Node Classes", fontsize=14)
for tick in ax2.xaxis.get_major_ticks():
    tick.label.set_fontsize(13)


add_stat_annotation(ax2, **params, data=dfMe[dfMe["Condition"] == "DHT"], box_pairs=boxPairs, test='Mann-Whitney')


fig.savefig(f"{figureRoot}/H3K27ac.pdf")
