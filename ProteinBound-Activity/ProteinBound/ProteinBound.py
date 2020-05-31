





# tfsLncap = readBed(home + "ALL.Prs.50.AllAg.LNCAP.bed")
# tfsLncap = sortBed(tfsLncap)

home2 = "/kuacc/users/ualtintas20"

conBed = readBed(home2 + "/regions/cons-arbs.bed")
conBed = sortBed(conBed)


indBed = readBed(home2 + "/regions/ind-arbs.bed")
indBed = sortBed(indBed)



nonBed = readBed(home2 + "/regions/Non-Active-ARBS.bed")
nonBed = sortBed(nonBed)


# ProteinBeds = {"arp":arpBed,
#                "ctc":ctcBed,
#                "fox":foxBed,
#                "med":medBed,
#                "ezh":ezhBed,
#                "ari":ariBed,
#                "pol":polBed,
#                "nmy":nmyBed
#
# }






home = "/kuacc/users/ualtintas20/LNCaP/BED"


proBed = readBed(home + "/TSS.hg19.Idx.bed")
proBed = sortBed(proBed)

# tsPBed = readBed(home + "/tsPActivity.bed")
# tsPBed = sortBed(tsPBed)
#
#
# tsMBed = readBed(home + "/tsMActivity.bed")
# tsMBed = sortBed(tsMBed)


# Histones = os.listdir(home + "/Histones")
# HistoneBeds = {}
# for histone in Histones:
#     histoneBed = readBed(home + "/Histones/" + histone)
#     histone = histone.split(".")[0]
#     histoneBed = sortBed(histoneBed)
#     HistoneBeds[histone] = histoneBed


Proteins = os.listdir(home + "/Proteins")
ProteinBeds = {}
for protein in Proteins:
    proteinBed = readBed(home + "/Proteins/" + protein)
    protein = protein.split(".")[0]
    proteinBed = sortBed(proteinBed)
    ProteinBeds[protein] = proteinBed



mappingValCon, mappingsMatrixCon = oneWayRangeSearch(conBed, ProteinBeds)

mappingValInd, mappingsMatrixInd = oneWayRangeSearch(indBed, ProteinBeds)

mappingValNon, mappingsMatrixNon = oneWayRangeSearch(nonBed, ProteinBeds)

mappingValPro, mappingsMatrixPro = oneWayRangeSearch(proBed, ProteinBeds)






conDf = pd.DataFrame.from_dict(mappingsMatrixCon, orient="index", columns=list(ProteinBeds.keys()))
conDf["NodeClass"] = "con"

indDf = pd.DataFrame.from_dict(mappingsMatrixInd, orient="index", columns=list(ProteinBeds.keys()))
indDf["NodeClass"] = "ind"

nonDf = pd.DataFrame.from_dict(mappingsMatrixNon, orient="index", columns=list(ProteinBeds.keys()))
nonDf["NodeClass"] = "non"

proDf = pd.DataFrame.from_dict(mappingsMatrixPro, orient="index", columns=list(ProteinBeds.keys()))
proDf["NodeClass"] = "pro"


BoundDf = pd.concat([conDf, indDf, nonDf, proDf])

pickle.dump(BoundDf, open("BoundDf.p", "wb"))


BoundDf = pickle.load(open("BoundDf.p", "rb"))
boundDf = BoundDf.drop(columns=["NodeClass"])
ArbsDf = BoundDf[BoundDf["NodeClass"] != "pro"]
arbsDf = ArbsDf.drop(columns=["NodeClass"])


lut = dict(zip(ArbsDf["NodeClass"].unique(), ["red","pink", "green", "yellow"]))
row_colors = ArbsDf["NodeClass"].map(lut)


fig = plt.figure(figsize=[20,50])
sns.clustermap(arbsDf, row_colors=row_colors, col_cluster=False)
plt.show()

fig = plt.figure(figsize=[20,50])
sns.heatmap(boundDf,cmap="vlag")
plt.savefig("b.pdf")



# mappingsProteinsCon = mappingsProteins
# nodeClass = ["con", "ind", "non"]
# data = dict.fromkeys(nodeClass, {})

dataCon = {}
for cre in mappingsMatrixCon.keys():
    dataCon[cre] = sum(mappingsMatrixCon[cre])

dataInd = {}
for cre in mappingsMatrixInd.keys():
    dataInd[cre] = sum(mappingsMatrixInd[cre])

dataNon = {}
for cre in mappingsMatrixNon.keys():
    dataNon[cre] = sum(mappingsMatrixNon[cre])


dataPro = {}
for cre in mappingsMatrixPro.keys():
    dataPro[cre] = sum(mappingsMatrixPro[cre])




data = {
    "non": list(dataNon.values()),
    "ind": list(dataInd.values()),
    "con": list(dataCon.values()),
    "pro": list(dataPro.values())
}


pickle.dump(data, open("mappingData.p", "wb"))


data = pickle.load(open("mappingData.p", "rb"))

combs = [[0, 2], [0, 1], [1, 2]]




combs = [
[0, 3],
[0, 2], [0, 1],
[1, 3], [1, 2], [2, 3]]

yDec = 1
y1 = 19
ylim = [0,20]


colorPalette = ["#d4f3ef", "#abf0e9", "#63b7af", "#ee8572"]

fig = plt.figure(figsize=(6,9))
boxPlot(data=data, combs=combs, colorPalette=colorPalette, y1=y1, yDec=yDec, ylim=ylim, paired=False)
plt.ylabel("# of Bound Protein")

fig.savefig("BoundProtein.pdf")
