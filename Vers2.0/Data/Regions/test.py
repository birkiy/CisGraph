import csv


def readBed(bed, file):
    with open(file) as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        for row in reader:
            # row3 = row[3].split(".")[0]
            row3 = row[3]
            bed[row3] = (row[0], int(row[1]), int(row[2]))

import bisect as bi
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


def sortedMapping(creBed, Beds):
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

            nTuple = rangesFromUpperRange(rangeX[0], rangeX[1], rangeX[2], Beds=beds)
            if len(nTuple) != 0:
                print(nTuple, rangeX)
                mappings[cre] = nTuple
            # if len(nTuple) >= 3 :
                plt.plot([rangeX[1], rangeX[2]], [0, 0], linewidth=2, color='k')
                plt.text((float(rangeX[1]) + float(rangeX[2])) * .5, 0.1, "cre", ha='center', va='bottom', color="k", fontsize=10)
                plt.xlim([rangeX[1]-10000, rangeX[2]+10000])
                plt.ylim([-0.5, len(nTuple)+0.5])
                y=1
                for tup in nTuple:
                    rangex = Beds[tup[1]][tup[0]]
                    plt.plot([rangex[1], rangex[2]], [y, y], linewidth=2, color='k')
                    plt.xlim([rangeX[1]-10000, rangeX[2]+10000])
                    plt.text((float(rangex[1]) + float(rangex[2])) * .5, y+0.1, tup[1], ha='center', va='bottom', color="k", fontsize=10)
                    y+=1
                plt.savefig("genes/" + cre + ".pdf")
                plt.close("all")
    return mappings



def sortBed(Gbed):
    Gbed = {k: v for k,v in sorted(Gbed.items(), key=lambda kv: kv[1][1])}
    Gbed = {k: v for k, v in sorted(Gbed.items(), key=lambda kv: kv[1][0])}
    return Gbed


creBed = {}
readBed(creBed, "creDHS.bed")
creBed = sortBed(creBed)
cre1 = {key: value for key, value in creBed.items() if value[0] == "chr1"}

proBed = {}
readBed(proBed, "promoters_ann_5kb.bed")
proBed = sortBed(proBed)
pro1 = {key: value for key, value in proBed.items() if value[0] == "chr1"}



conBed = {}
readBed(conBed, "cons-arbs.bed")
conBed = sortBed(conBed)
con1 = {key: value for key, value in conBed.items() if value[0] == "chr1"}


indBed = {}
readBed(indBed, "ind-arbs.bed")
indBed = sortBed(indBed)
ind1 = {key: value for key, value in indBed.items() if value[0] == "chr1"}


nonBed = {}
readBed(nonBed, "Non-Active-ARBS.bed")
nonBed = sortBed(nonBed)
non1 = {key: value for key, value in nonBed.items() if value[0] == "chr1"}



arpBed = {}
readBed(arpBed, "ProteinBounds/AR.bed")
sortBed(arpBed)
arpBed = {k: v for k,v in sorted(arpBed.items(), key=lambda kv: kv[1][1])}

ctcBed = {}
readBed(ctcBed, "ProteinBounds/CTCF.bed")
sortBed(ctcBed)
ctcBed = {k: v for k,v in sorted(ctcBed.items(), key=lambda kv: kv[1][1])}

foxBed = {}
readBed(foxBed, "ProteinBounds/FOXA1.bed")
sortBed(foxBed)
foxBed = {k: v for k,v in sorted(foxBed.items(), key=lambda kv: kv[1][1])}

medBed = {}
readBed(medBed, "ProteinBounds/MED1.bed")
sortBed(foxBed)
medBed = {k: v for k,v in sorted(medBed.items(), key=lambda kv: kv[1][1])}

enzBed = {}
readBed(enzBed, "ProteinBounds/EZH2.bed")
sortBed(enzBed)
enzBed = {k: v for k,v in sorted(enzBed.items(), key=lambda kv: kv[1][1])}

ariBed = {}
readBed(ariBed, "ProteinBounds/ARID.bed")
sortBed(ariBed)
ariBed = {k: v for k,v in sorted(ariBed.items(), key=lambda kv: kv[1][1])}

polBed = {}
readBed(polBed, "ProteinBounds/POLR2.bed")
sortBed(polBed)
polBed = {k: v for k,v in sorted(polBed.items(), key=lambda kv: kv[1][1])}


nmyBed = {}
readBed(nmyBed, "ProteinBounds/NMYC.bed")
sortBed(nmyBed)
nmyBed = {k: v for k,v in sorted(nmyBed.items(), key=lambda kv: kv[1][1])}




Beds = {"arp":arpBed,
        "ctc":ctcBed,
        "fox":foxBed,
        "med":medBed,
        "enz":enzBed,
        "ari":ariBed,
        "pol":polBed,
        "nmy":nmyBed,

        "pro": proBed,
        "con": conBed,
        "ind": indBed,
        "non": nonBed

}

#
# Beds = {"pro": proBed,
#         "con": conBed,
#         "ind": indBed,
#         "non": nonBed}
import matplotlib.pyplot as plt

import time
start_time = time.time()
mappings = sortedMapping(creBed, Beds)

# mappings = sortedMapping(proBed, Beds)
print("--- %s seconds ---" % (time.time() - start_time))

import pickle
#mappings = pickle.load(open("mappings.p", "rb"))

ones = {k: v for k, v in mappings.items() if len(v) == 1}
twos = {k: v for k, v in mappings.items() if len(v) == 2}
threes = {k: v for k, v in mappings.items() if len(v) == 3}

fours = {k: v for k, v in mappings.items() if len(v) == 4}



print(len(ones.keys()))
print(len(twos.keys()))
print(len(threes.keys()))
print(len(fours.keys()))

print(fours)

pickle.dump(mappings, open("mappings.p", "wb"))


#
#
# mapping = {}
# for cre in creBed:
#     rangeX = creBed[cre]
#
#     nTuple = rangesFromUpperRange(rangeX[0], rangeX[1], rangeX[2], Beds=Beds)
#     if len(nTuple) != 0:
#         print(nTuple, rangeX)
#         mapping[cre] = nTuple
#     if len(nTuple) == 3 :
#         plt.plot([rangeX[1], rangeX[2]], [0, 0], linewidth=2, color='k')
#         plt.text((float(rangeX[1]) + float(rangeX[2])) * .5, 0.1, "cre", ha='center', va='bottom', color="k", fontsize=10)
#         plt.xlim([rangeX[1]-10000, rangeX[2]+10000])
#         y=1
#         for tup in nTuple:
#             rangex = Beds[tup[1]][tup[0]]
#             plt.plot([rangex[1], rangex[2]], [y, y], linewidth=2, color='k')
#             plt.text((float(rangex[1]) + float(rangex[2])) * .5, y+0.1, tup[1], ha='center', va='bottom', color="k", fontsize=10)
#             y+=1
#         plt.savefig(cre + ".pdf")
#
#
#
#
#


#
# #
#
#
# #
# #
# # # "cre": cre1,"pro": pro1,
# # chr1 = {
# #
# #         "con": con1,
# #         "ind": ind1,
# #         "non": non1
# # }
# #
# #
# #
# # y = 0
# # for c1 in chr1:
# #     for reg in chr1[c1]:
# #         print(reg, y)
# #         x1, x2 = chr1[c1][reg][1], chr1[c1][reg][2]
# #         plt.plot([x1, x2], [y, y], linewidth=1, color='k')
# #
# #
# #     y += 1
# #
# # plt.savefig("test.pdf")
