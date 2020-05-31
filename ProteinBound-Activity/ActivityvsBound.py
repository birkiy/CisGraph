


home = "/kuacc/users/ualtintas20/LNCaP/BED"


tsPBed = readBed(home + "/tsPActivity.bed")
tsPBed = sortBed(tsPBed)


tsMBed = readBed(home + "/tsMActivity.bed")
tsMBed = sortBed(tsMBed)





arbsBed = readBed(home + "/arbsActivity.bed")
arbsBed = sortBed(arbsBed)



Proteins = os.listdir(home + "/Proteins")
ProteinBeds = {}
for protein in Proteins:
    proteinBed = readBed(home + "/Proteins/" + protein)
    protein = protein.split(".")[0]
    proteinBed = sortBed(proteinBed)
    ProteinBeds[protein] = proteinBed


mappingValTsP, mappingsMatrixTsP = oneWayRangeSearch(tsPBed, ProteinBeds)

mappingValTsM, mappingsMatrixTsM = oneWayRangeSearch(tsMBed, ProteinBeds)



mappingValArbs, mappingsMatrixArbs = oneWayRangeSearch(arbsBed, ProteinBeds)


names = ["#chrom","start","end",	"name",	"score",	"strand",	"thickStart",	"thickEnd",	"itemRGB",	"blockCount",	"blockSizes",	"blockStart", "deepTools_group","groseq_dht.plus","groseq_dht.plus.1","groseq_dht.minus","groseq_dht.minus.1","groseq_dmso.plus","groseq_dmso.plus.1","groseq_dmso.minus","groseq_dmso.minus.1", "Strand", "nodeClass", "Activity"]
arbsActivity = pd.read_csv(home + "/arbsActivity.bed", sep="\t", names=names)

BoundTFarbs = {}
for key in mappingsMatrixArbs.keys():
    cre = key[0]
    BoundTFarbs[cre] = mappingsMatrixArbs[key]

dataArbs = []
for cre in list(arbsActivity["name"]):
    dataArbs += [sum(BoundTFarbs[cre])]

arbsActivity["BoundTF"] = dataArbs

arbsActivity.to_csv("ArbsActivityvsBound.tab", sep="\t", index=False, header=False)




names = ["#chrom","start","end",	"name",	"score",	"strand",	"thickStart",	"thickEnd",	"itemRGB",	"blockCount",	"blockSizes",	"blockStart", "deepTools_group","groseq_dht.plus","groseq_dht.plus.1","groseq_dht.minus","groseq_dht.minus.1","groseq_dmso.plus","groseq_dmso.plus.1","groseq_dmso.minus","groseq_dmso.minus.1", "Strand", "nodeClass", "Activity"]

arbsActivity = pd.read_csv("ArbsActivityvsBound.tab", sep="\t", names=names + ["BoundTF"])

arbsActivity.sort_values(by=["BoundTF"])
arbsActivity.reset_index()
# arbsActivity["quantile"] = pd.qcut(arbsActivity["BoundTF"], q=)



y = "Activity"
x = "nodeClass"
# y = "BoundTF"
hue = "Strand"
hue_order=['+', '-']

fig = plt.figure(figsize=[9,10])
colorPalette = ["#63b7af", "#abf0e9", "#d4f3ef"]
colorPalette = ["#FFA96B", "#75E3C1"]
ax = sns.boxplot(y=y, x=x, hue=hue, data=arbsActivity, palette=colorPalette)
ax.set_ylim([0,10])



box_pairs = box_pairs + [((clar, '+'), (clar, '-')) for clar in arbsActivity['nodeClass'].unique()]

box_pairs=[
    (("con", "+"), ("con", "-")),
    (("ind", "+"), ("ind", "-")),
    (("non", "+"), ("non", "-")),
    (("con", "+"), ("ind", "+")),
    (("con", "+"), ("non", "+")),
    (("ind", "+"), ("non", "+")),
    (("con", "-"), ("ind", "-")),
    (("con", "-"), ("non", "-")),
    (("ind", "-"), ("non", "-")),
    ]

add_stat_annotation(ax, data=arbsActivity, x=x, y=y, hue=hue, loc="inside", test="Mann-Whitney", box_pairs=box_pairs)

plt.savefig('example_barplot_hue.png', dpi=300, bbox_inches='tight')


combs = [[0, 3],
[1, 3],[2, 3], [0, 2],
[1, 2], [0, 1]]

PVal =[]
for c in combs:
    if paired:
        pVal = scipy.stats.wilcoxon(data[list(data.keys())[c[0]]], data[list(data.keys())[c[1]]])[1]
    else:
        pVal = scipy.stats.ranksums(data[list(data.keys())[c[0]]], data[list(data.keys())[c[1]]])[1]
    PVal += [pVal]
if multiTest:
    pAdj = multipletests(PVal, method='bonferroni')
else:
    pAdj = PVal


plt.ylabel("ARBSs' Promoter Activty")
plt.xlabel("# of Bound Protein")
plt.savefig("ArbsBoundvsActivity.pdf")


conActivity = arbsActivity[arbsActivity["nodeClass"] == "con"]
scipy.stats.spearmanr(conActivity["Activity"], conActivity["BoundTF"])

colorPalette = ["#FFA96B", "#75E3C1"]
sns.boxplot(y="Activity", x="nodeClass", hue="Strand", data=conActivity, palette=colorPalette)













tsP = pd.read_csv(home + "/tsPActivity.bed", sep="\t", names=["#chrom", "start", "end", "name", "Activity"])

BoundTFtsP = {}
for key in mappingsMatrixTsP.keys():
    cre = key[0]
    BoundTFtsP[cre] = mappingsMatrixTsP[key]

dataTsP = []
for cre in list(tsP["name"]):
    dataTsP += [sum(BoundTFtsP[cre])]

tsP["BoundTF"] = dataTsP
tsP["Strand"] = "+"

tsP.to_csv("tsPActivityvsBound.tab", sep="\t", index=False, header=False)







tsM = pd.read_csv(home + "/tsMActivity.bed", sep="\t", names=["#chrom", "start", "end", "name", "Activity"])

BoundTFtsM = {}
for key in mappingsMatrixTsM.keys():
    cre = key[0]
    BoundTFtsM[cre] = mappingsMatrixTsM[key]

dataTsM = []
for cre in list(tsM["name"]):
    dataTsM += [sum(BoundTFtsM[cre])]

tsM["BoundTF"] = dataTsM
tsM["Strand"] = "-"
tsM.to_csv("tsMActivityvsBound.tab", sep="\t", index=False, header=False)


home = "/home/birkiy/github/CisGraph/GRO"


tsM = pd.read_csv(home + "/tsMActivityvsBound.tab", sep="\t", names=["#chrom", "start", "end", "name", "Activity", "BoundTF"])
tsM["Strand"] = "-"
tsP = pd.read_csv(home + "/tsPActivityvsBound.tab", sep="\t", names=["#chrom", "start", "end", "name", "Activity", "BoundTF"])
tsP["Strand"] = "+"
tss = pd.concat([tsP, tsM])


home = "/home/birkiy/github/CisGraph/GRO"
# tss = pd.read_csv(home + "/tssActivityvsBound.tab", sep="\t", names=["#chrom", "start", "end", "name", "Activity"])


tss.sort_values(by=["BoundTF"])
tss.reset_index()
tss["quantile"] = pd.qcut(tss["BoundTF"], q=)

fig = plt.figure(figsize=[9,10])
colorPalette = ["#FFA96B", "#75E3C1"]
sns.boxplot(y="Activity", x="BoundTF", hue="Strand", data=tss, palette=colorPalette)
plt.ylabel("Promoter Activty Quantiles")
plt.xlabel("# of Bound Protein")
plt.savefig("BoundvsActivity.pdf")

n = len(tss)




range(0,23, 2)




scipy.stats.wilcoxon(tss["Activity"], tss["BoundTF"])

ax =
