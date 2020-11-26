


from Functions.Packages import *

chip = pd.read_csv(f"{dataRoot}/ARBS_chip_seq_2.csv")




conAnn = pd.read_csv(f"{dataRoot}/Regions/con.AnnFile.cvs", sep="\t")
indAnn = pd.read_csv(f"{dataRoot}/Regions/ind.AnnFile.csv", sep="\t")
nonAnn = pd.read_csv(f"{dataRoot}/Regions/non.AnnFile.csv", sep="\t")

conAnn = conAnn.rename(columns={"seqnames": "Chr", "start": "Start", "end": "End"})
indAnn = indAnn.rename(columns={"seqnames": "Chr", "start": "Start", "end": "End"})
nonAnn = nonAnn.rename(columns={"seqnames": "Chr", "start": "Start", "end": "End"})




con = conAnn[["V4", "Chr", "Start", "End", "annotation", "geneStrand", "SYMBOL"]]
ind = indAnn[["V4", "Chr", "Start", "End", "annotation", "geneStrand", "SYMBOL"]]
non = nonAnn[["V4", "Chr", "Start", "End", "annotation", "geneStrand", "SYMBOL"]]


con["nodeClass"] = "con"
ind["nodeClass"] = "ind"
non["nodeClass"] = "non"


conChip = chip[chip["Unnamed: 0"].isin(con["V4"])]
indChip = chip[chip["Unnamed: 0"].isin(ind["V4"])]
nonChip = chip[chip["Unnamed: 0"].isin(non["V4"])]

chip["nodeClass"] = "NA"

chip.loc[chip["Unnamed: 0"].isin(con["V4"]), "nodeClass"] = "con"
chip.loc[chip["Unnamed: 0"].isin(ind["V4"]), "nodeClass"] = "ind"
chip.loc[chip["Unnamed: 0"].isin(non["V4"]), "nodeClass"] = "non"



hMeMark = chip[["Unnamed: 0", "nodeClass", "H3K4me3-ETOH-GSE114732", "H3K4me3-R1881-GSE114732", "H3K4me1-ETOH-GSE114732", "H3K4me1-R1881-GSE114732"]]

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



#
#
# fig = plt.figure(figsize=[8,9])
# plt.subplots_adjust(top = 0.75)
#
# params = dict(data=dfMe,
#               x='nodeClass',
#               y="H3K4me3/H3K4me1",
#               hue='Condition',
#               order=["con", "ind", "non", "NA"],
#               hue_order=["EtOH", "R1881"])
#
#
#
# ax1 = sns.boxplot( **params)
# ax1.set(yscale="log")
#
# boxPairs = [
#     (("con", "EtOH"), ("con", "R1881")),
#     (("ind", "EtOH"), ("ind", "R1881")),
#     (("non", "EtOH"), ("non", "R1881")),
#     (("NA", "EtOH"), ("NA", "R1881"))
#
# ]
# add_stat_annotation(ax1, **params, box_pairs=boxPairs, test='Mann-Whitney', loc='inside')
# boxPairs = [
#     (("con", "EtOH"), ("ind", "EtOH")),
#     (("con", "EtOH"), ("non", "EtOH")),
#     (("con", "EtOH"), ("NA", "EtOH")),
#     (("ind", "EtOH"), ("non", "EtOH")),
#     (("ind", "EtOH"), ("NA", "EtOH")),
#     (("non", "EtOH"), ("NA", "EtOH")),
#
#     (("con", "R1881"), ("ind", "R1881")),
#     (("con", "R1881"), ("non", "R1881")),
#     (("con", "R1881"), ("NA", "R1881")),
#     (("ind", "R1881"), ("non", "R1881")),
#     (("ind", "R1881"), ("NA", "R1881")),
#     (("non", "R1881"), ("NA", "R1881"))
#
# ]
# add_stat_annotation(ax1, **params, box_pairs=boxPairs, test='Mann-Whitney', loc='outside')
#
#
#
# fig.savefig(f"{figureRoot}/H3K4me3-H3K4me1.ARBS.pdf")
#
#
#


fig = plt.figure(figsize=[9, 12])
gs = gridspec.GridSpec(ncols=2, nrows=2)
# plt.subplots_adjust(top = 0.8, hspace=0.5)

colorPalette = ["#63b7af", "#abf0e9","#d4f3ef", "#f5fffd", "#ee8572"]

params = dict(x='nodeClass',
              order=["con", "ind", "non", "NA"])

boxPairs = [("NA", "con"),
            ("NA", "ind"),
            ("NA", "non"),
            ("con", "ind"),
            ("con", "non"),
            ("ind", "non")
            ]

ax1 = fig.add_subplot(gs[0])
ax1 = sns.boxplot(data=dfMe[dfMe["Condition"] == "EtOH"], y="H3K4me3-GSE114732", palette=colorPalette, **params)
plt.title("EtOH")
add_stat_annotation(ax1, data=dfMe[dfMe["Condition"] == "EtOH"], y="H3K4me3-GSE114732", **params, box_pairs=boxPairs, test='Mann-Whitney')

ax2 = fig.add_subplot(gs[1])
ax2 = sns.boxplot(data=dfMe[dfMe["Condition"] == "R1881"], y="H3K4me3-GSE114732", palette=colorPalette, **params)
plt.title("R1881")
add_stat_annotation(ax2,data=dfMe[dfMe["Condition"] == "R1881"], y="H3K4me3-GSE114732", **params, box_pairs=boxPairs, test='Mann-Whitney')


ax3 = fig.add_subplot(gs[2])
ax3 = sns.boxplot(data=dfMe[dfMe["Condition"] == "EtOH"], y="H3K4me1-GSE114732", palette=colorPalette, **params)
plt.title("EtOH")
add_stat_annotation(ax3, data=dfMe[dfMe["Condition"] == "EtOH"], y="H3K4me1-GSE114732", **params, box_pairs=boxPairs, test='Mann-Whitney')

ax4 = fig.add_subplot(gs[3])
ax4 = sns.boxplot(data=dfMe[dfMe["Condition"] == "R1881"], y="H3K4me1-GSE114732", palette=colorPalette, **params)
plt.title("R1881")
add_stat_annotation(ax4,data=dfMe[dfMe["Condition"] == "R1881"], y="H3K4me1-GSE114732", **params, box_pairs=boxPairs, test='Mann-Whitney')


fig.savefig(f"{figureRoot}/H3K4me3.H3K4me1.ARBS.pdf")





rnaPol = chip[["Unnamed: 0", "nodeClass", 'POLR2A-ETOH-GSE43253', "POLR2A-DHT-GSE43253"]]


etMe = rnaPol[["Unnamed: 0", "nodeClass", 'POLR2A-ETOH-GSE43253']]
etMe["Condition"] = "EtOH"
etMe = etMe.rename(columns={'POLR2A-ETOH-GSE43253': 'POLR2A-GSE43253'})



r1Me = rnaPol[["Unnamed: 0", "nodeClass", "POLR2A-DHT-GSE43253"]]
r1Me["Condition"] = "DHT"
r1Me = r1Me.rename(columns={'POLR2A-DHT-GSE43253': 'POLR2A-GSE43253'})


dfMe = pd.concat([etMe, r1Me])



fig = plt.figure(figsize=[9, 6])
gs = gridspec.GridSpec(ncols=2, nrows=1)
# plt.subplots_adjust(top = 0.8, hspace=0.5)


colorPalette = ["#63b7af", "#abf0e9","#d4f3ef", "#f5fffd", "#ee8572"]

params = dict(x='nodeClass',
              y='POLR2A-GSE43253',
              order=["con", "ind", "non", "NA"])

boxPairs = [("NA", "con"),
            ("NA", "ind"),
            ("NA", "non"),
            ("con", "ind"),
            ("con", "non"),
            ("ind", "non")
            ]

ax1 = fig.add_subplot(gs[0])
ax1 = sns.boxplot(**params, data=dfMe[dfMe["Condition"] == "EtOH"], palette=colorPalette)
plt.title("EtOH")
add_stat_annotation(ax1, **params, data=dfMe[dfMe["Condition"] == "EtOH"], box_pairs=boxPairs, test='Mann-Whitney')

ax2 = fig.add_subplot(gs[1])
ax2 = sns.boxplot(**params, data=dfMe[dfMe["Condition"] == "DHT"], palette=colorPalette)
plt.title("DHT")
add_stat_annotation(ax2, **params, data=dfMe[dfMe["Condition"] == "DHT"], box_pairs=boxPairs, test='Mann-Whitney')


#
# boxPairs = [
#     (("con", "EtOH"), ("ind", "EtOH")),
#     (("con", "EtOH"), ("non", "EtOH")),
#     (("con", "EtOH"), ("NA", "EtOH")),
#     (("ind", "EtOH"), ("non", "EtOH")),
#     (("ind", "EtOH"), ("NA", "EtOH")),
#     (("non", "EtOH"), ("NA", "EtOH")),
#
#     (("con", "DHT"), ("ind", "DHT")),
#     (("con", "DHT"), ("non", "DHT")),
#     (("con", "DHT"), ("NA", "DHT")),
#     (("ind", "DHT"), ("non", "DHT")),
#     (("ind", "DHT"), ("NA", "DHT")),
#     (("non", "DHT"), ("NA", "DHT"))
#
# ]
# boxPairs = [((x[0][1], x[0][0]), (x[1][1], x[1][0])) for x in boxPairs]
#
# add_stat_annotation(ax1, **params, box_pairs=boxPairs, test='Mann-Whitney', loc='inside')
#
#
# boxPairs = [
#     (("con", "EtOH"), ("con", "DHT")),
#     (("ind", "EtOH"), ("ind", "DHT")),
#     (("non", "EtOH"), ("non", "DHT")),
#     (("NA", "EtOH"), ("NA", "DHT"))
#
# ]
# boxPairs = [((x[0][1], x[0][0]), (x[1][1], x[1][0])) for x in boxPairs]
#
# add_stat_annotation(ax1, **params, box_pairs=boxPairs, test='Mann-Whitney', loc='outside')


fig.savefig(f"{figureRoot}/RNAPII.ARBS.pdf")
