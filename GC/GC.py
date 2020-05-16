

import csv
import pickle
import bisect as bi
import matplotlib.pyplot as plt
import time
import scipy
import pandas as pd
import seaborn as sns



def rangesFromUpperRange(chrKey, start, end, Beds):
    values = []
    beds = {}
    for Bed in Beds:
        bed = {key: value for key, value in Beds[Bed].items() if value[0] == chrKey}
        bed = sortBed(bed)
        beds[Bed] = bed
    bedName = list(beds.keys())
    beds = list(beds.values())
    for idx, enhancer in enumerate(beds):
        enhancerEnd = [(value[2], value[1], value[0], key) for key, value in enhancer.items()]
        left = bi.bisect_left(enhancerEnd, (start, end))

        enhancerStart = [(value[1], value[2], value[0], key) for key, value in enhancer.items()]
        right = bi.bisect_right(enhancerStart, (end, start))
        values += [
            (keys, bedName[idx])
            for keys, value in list(enhancer.items())[left:right]
            if chrKey == value[0] and
               value[2] - start > start - value[1] and
               end - value[1] > value[2] - end
        ]
    return values


def sortedMapping(creBed, Beds, pro=False):
    chrs = set([cre[0] for cre in creBed.values()])

    mappings = {}
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
            if pro:
                minus = 3000
            else:
                minus = 0

            nTuple = rangesFromUpperRange(rangeX[0], rangeX[1]+minus, rangeX[2]-minus, Beds=beds)
            if len(nTuple) != 0:
                # print(nTuple, rangeX)
                # mappings[cre] = nTuple

                C = [GCmetBed[met][3] for met, i in nTuple]
                cov = [GCmetBed[met][4] for met, i in nTuple]

                mappings[cre] = (sum(C) / sum(cov))
                # mappings[cre] = len(cov)
            #
            #
            # # if len(nTuple) >= 3 :
            #     plt.plot([rangeX[1], rangeX[2]], [0, 0], linewidth=2, color='k')
            #     plt.text((float(rangeX[1]) + float(rangeX[2])) * .5, 0.1, "ind", ha='center', va='bottom', color="k", fontsize=10)
            #     plt.xlim([rangeX[1]-1000, rangeX[2]+1000])
            #     plt.ylim([-0.5, len(nTuple)+0.5])
            #     y=1
            #     for tup in nTuple:
            #         rangex = Beds[tup[1]][tup[0]]
            #         plt.plot([rangex[1], rangex[2]], [y, y], linewidth=2, color='k')
            #         plt.xlim([rangeX[1]-1000, rangeX[2]+1000])
            #         plt.text((float(rangex[1]) + float(rangex[2])) * .5, y+0.1, tup[1], ha='center', va='bottom', color="k", fontsize=10)
            #         y+=1
            #     plt.savefig("ind/" + cre + ".pdf")
            #     plt.close("all")
    return mappings

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


            # if not geneSymbol in geneNodes:
            #     continue

            if (logFC == "NA" or Qval == "NA") or (not geneSymbol in bed.keys()):
                continue

            logFC = float(logFC)
            Qval = float(Qval)


            if (logFC > 1 and Qval < 0.05):
                upP[geneSymbol] = bed[geneSymbol]

            elif (logFC < -1 and Qval < 0.05):
                dwP[geneSymbol] = bed[geneSymbol]
    return upP, dwP



def readBed(bed, file):
    with open(file) as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        for row in reader:
            # row3 = row[3].split(".")[0]
            row3 = row[3]
            bed[row3] = (row[0], int(row[1]), int(row[2]))

def sortBed(Gbed):
    Gbed = {k: v for k,v in sorted(Gbed.items(), key=lambda kv: kv[1][1])}
    Gbed = {k: v for k, v in sorted(Gbed.items(), key=lambda kv: kv[1][0])}
    return Gbed


def flat(list):
    return [val for sub in list for val in sub]


proBed = {}
readBed(proBed,"TSS.hg19.+.bed")
proBed = sortBed(proBed)

upPBed, dwPBed = logFC("GSE64529_diffexpr-results.csv", proBed)
f= open("TSS.hg19.+.up.bed","w+")
for region in upPBed:
    f.write(upPBed[region][0]  + "\t" + str(upPBed[region][1]) + "\t" + str(upPBed[region][2]) + "\t" + region + "\n")
f.close()

f= open("TSS.hg19.+.dw.bed","w+")
for region in dwPBed:
    f.write(dwPBed[region][0]  + "\t" + str(dwPBed[region][1]) + "\t" + str(dwPBed[region][2]) + "\t" + region + "\n")
f.close()


# conBed = {}
# readBed(conBed,"cons-arbs.bed")
# conBed = sortBed(conBed)
#
# indBed = {}
# readBed(indBed,"ind-arbs.bed")
# indBed = sortBed(indBed)
#
#
# nonBed = {}
# readBed(nonBed, "Non-Active-ARBS.bed")
# nonBed = sortBed(nonBed)
#
#
# GCmetBed = {}
# with open("GCint.bed") as tsvfile:
#     reader = csv.reader(tsvfile, delimiter='\t')
#     for row in reader:
#         # row3 = row[3].split(".")[0]
#         row3 = row[3]
#         GCmetBed[row3] = (row[0], int(row[1]), int(row[2]), int(row[4]), int(row[5]))
# GCmetBed = sortBed(GCmetBed)
# Beds = {"met": GCmetBed}
#
# start_time = time.time()
# # mappingsDwP = sortedMapping(dwPBed, Beds)
# # mappingsUpP = sortedMapping(upPBed, Beds)
# mappingsPro = sortedMapping(proBed, Beds, pro=True)
# mappingsCon = sortedMapping(conBed, Beds)
# mappingsInd = sortedMapping(indBed, Beds)
# mappingsNon = sortedMapping(nonBed, Beds)
# #"#27496d",
# # mappings = {"dwP": list(mappingsDwP.values()),
# #             "upP": list(mappingsUpP.values()),
# #             "con": list(mappingsCon.values()),
# #             "ind": list(mappingsInd.values()),
# #             "non": list(mappingsNon.values())}
# #
#
#
# mappings = {"pro": list(mappingsPro.values()),
#             "con": list(mappingsCon.values()),
#             "ind": list(mappingsInd.values()),
#             "non": list(mappingsNon.values())}
#
#
#
# cGC = flat([mappings[d] for d in mappings])
# nodeClasses = [nodeClass for nodeClass in mappings.keys() for _ in range(len(mappings[nodeClass]))]
#
#
# df = {
#     "nodeClass": nodeClasses,
#     "GC": cGC
# }
#
# df = pd.DataFrame(df)
#
# data = mappings
# print(df)
# fig = plt.figure(figsize=(6,9))
# sns.boxplot(x="nodeClass", y="GC", data=df, palette=[ "#ee8572", "#63b7af","#d4f3ef", "#abf0e9"])
# #
# # combs = [[0, 4],
# # [0, 3], [1, 4],
# # [0, 2], [1, 3], [2, 4],
# # [0, 1], [1, 2], [2, 3], [3, 4]]
#
#
#
# combs = [
# [0, 3],
# [0, 2], [1, 3],
# [0, 1], [1, 2], [2, 3]]
# yDec = 0.05
# y1 = 2
# for c in combs:
#     pVal = scipy.stats.ranksums(data[list(data.keys())[c[0]]], data[list(data.keys())[c[1]]])[1]
#
#
#     x1 = c[0]
#     x2 = c[1]
#     y1 = y1 - yDec
#     y2 = y1 + (yDec / 3)
#
#
#     if pVal < 0.0001:
#         pS = "***"
#     elif pVal < 0.001:
#         pS = "**"
#     elif pVal < 0.01:
#         pS = "*"
#     elif pVal < 0.05:
#         pS = "."
#     else:
#         pS = "ns"
#
#     plt.plot([x1, x1, x2, x2], [y1, y2, y2, y1], linewidth=1, color='k')
#     plt.text((x1 + x2) * .5, y2, pS, ha='center', va='bottom', color="k", fontsize=10)
#
# plt.savefig("GCmetMean.pdf")
# plt.close("all")
#
# # print(mappings)
# # mappings = sortedMapping(proBed, Beds)
# print("--- %s seconds ---" % (time.time() - start_time))
#
#
#
# pickle.dump(GCmetBed, open("GCmet.p", "wb"))
# # print(GCmetBed)
