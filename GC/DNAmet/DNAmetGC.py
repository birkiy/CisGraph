
# Corrected tests





def sortedMapping(creBed, Beds):
    chrs = set([cre[0] for cre in creBed.values()])
    mappingsMean = {}
    mappingsLen = {}
    for chr1 in chrs:
        cre1 = {key: value for key, value in creBed.items() if value[0] == chr1}
        cre1 = sortBed(cre1)
        beds = dict.fromkeys(Beds.keys(), {})
        for bed in Beds:
            bed1 = {key: value for key, value in Beds[bed].items() if value[0] == chr1}
            bed1 = sortBed(bed1)
            beds[bed] = bed1
        for cre in cre1:
            rangeX = cre1[cre]
            nTuple = rangesFromUpperRange(rangeX[0], rangeX[1], rangeX[2], Beds=beds)
            if len(nTuple) != 0:
                C = [GCmetBed[met][3] for met, i in nTuple]
                cov = [GCmetBed[met][4] for met, i in nTuple]
                mappingsMean[cre] = (sum(C) / sum(cov))
                mappingsLen[cre] = len(cov)
    return mappingsMean, mappingsLen




def logFC(file, bed):
    upP = {}
    dwP = {}
    with open(file) as csvFile:
        reader = csv.reader(csvFile, delimiter=',')
        next(reader)
        for row in reader:
            geneSymbol = row[1]
            logFC = row[3]
            Qval = row[7]
            if (logFC == "NA" or Qval == "NA") or (not geneSymbol in bed.keys()):
                continue
            logFC = float(logFC)
            Qval = float(Qval)
            if (logFC > 1 and Qval < 0.05):
                upP[geneSymbol] = bed[geneSymbol]
            elif (logFC < -1 and Qval < 0.05):
                dwP[geneSymbol] = bed[geneSymbol]
    return upP, dwP



proBed = readBed("TSS.hg19.uniq.bed")

proBed = readBed("promoters_ann_5kb.bed")
proBed = sortBed(proBed)
#
# upPBed, dwPBed = logFC("GSE64529_diffexpr-results.csv", proBed)
# f= open("TSS.hg19.+.up.bed","w+")
# for region in upPBed:
#     f.write(upPBed[region][0]  + "\t" + str(upPBed[region][1]) + "\t" + str(upPBed[region][2]) + "\t" + region + "\n")
# f.close()
#
# f= open("TSS.hg19.+.dw.bed","w+")
# for region in dwPBed:
#     f.write(dwPBed[region][0]  + "\t" + str(dwPBed[region][1]) + "\t" + str(dwPBed[region][2]) + "\t" + region + "\n")
# f.close()

upPBed, dwPBed = logFC("GSE64529_diffexpr-results.csv", proBed)

conBed = readBed("cons-arbs.bed")
conBed = sortBed(conBed)


indBed = readBed("ind-arbs.bed")
indBed = sortBed(indBed)



nonBed = readBed("Non-Active-ARBS.bed")
nonBed = sortBed(nonBed)


GCmetBed = {}
with open("GCint.bed") as tsvfile:
    reader = csv.reader(tsvfile, delimiter='\t')
    for row in reader:
        # row3 = row[3].split(".")[0]
        row3 = row[3]
        GCmetBed[row3] = (row[0], int(row[1]), int(row[2]), int(row[4]), int(row[5]))


GCmetBed = sortBed(GCmetBed)
Beds = {"met": GCmetBed}



start_time = time.time()
# mappingsDwP = sortedMapping(dwPBed, Beds)
# mappingsUpP = sortedMapping(upPBed, Beds)
mappingsUpP = sortedMapping(upPBed, Beds)
mappingsDwP = sortedMapping(dwPBed, Beds)

mappingsProMean, mappingsProLen = sortedMapping(proBed, Beds)
mappingsConMean, mappingsConLen = sortedMapping(conBed, Beds)
mappingsIndMean, mappingsIndLen = sortedMapping(indBed, Beds)
mappingsNonMean, mappingsNonLen = sortedMapping(nonBed, Beds)
#"#27496d",
# mappings = {"dwP": list(mappingsDwP.values()),
#             "upP": list(mappingsUpP.values()),
#             "con": list(mappingsCon.values()),
#             "ind": list(mappingsInd.values()),
#             "non": list(mappingsNon.values())}
#



fig = plt.figure(figsize=[12,15])
gs = gridspec.GridSpec(ncols=2, nrows=2, width_ratios=[4, 4], height_ratios=[6, 3])

fig.add_subplot(gs[0,0])

colorPalette = ["#d4f3ef", "#abf0e9", "#63b7af","#ee8572"]

combs = [
[0, 3],
[0, 2], [0, 1],
[1, 3], [1, 2], [2, 3]]

yDec = 0.05
y1 = 1.3
ylim = [0,1.3]

data = {"non": list(mappingsNonMean.values()),
        "ind": list(mappingsIndMean.values()),
        "con": list(mappingsConMean.values()),
        "pro": list(mappingsProMean.values())}

boxPlot(data=data, combs=combs, colorPalette="Greys", y1=y1, yDec=yDec, ylim=ylim)


plt.ylabel("Average DNA methylation (C/Coverage)", fontsize=15)

fig.add_subplot(gs[0,1])

yDec = 1.5
y1 = 37
ylim = [0,38]

data = {"non": list(mappingsNonLen.values()),
        "ind": list(mappingsIndLen.values()),
        "con": list(mappingsConLen.values()),
        "pro": list(mappingsProLen.values())}

boxPlot(data=data, combs=combs, colorPalette="Greys", y1=y1, yDec=yDec, ylim=ylim)
plt.ylabel("# of DNA methylation per RE", fontsize=15)


fig.add_subplot(gs[1,:])
r = [0,1,2,3]

beds = {"non": list(nonBed.values()),
        "ind": list(indBed.values()),
        "con": list(conBed.values()),
        "pro": list(proBed.values())}


df = {
    "nodeClass": list(data.keys()),
    "Met": [len(v) for v in list(data.values())],
    "UnM": [len(b) - len(d) for b, d in zip( list(beds.values()), list(data.values()) )]
}
df = pd.DataFrame(df)

totals = [i+j for i,j in zip(df['Met'], df['UnM'])]
metC = [i / j * 100 for i,j in zip(df['Met'], totals)]
unmC = [i / j * 100 for i,j in zip(df['UnM'], totals)]
barWidth = 0.85
names = tuple(beds.keys())
plt.barh(r, metC, color='#5E5E5E', edgecolor='white')
# Create orange Bars
plt.barh(r, unmC ,left=metC, color='#A5A5A5', edgecolor='white')
plt.xlabel("% of Methylated Element", fontsize=15)
plt.legend(("Methylated", "Unmethylated"))
plt.yticks(r, names)

fig.savefig("DNAmet.pdf")
