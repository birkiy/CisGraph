
from MainFunctions.Utils import *


def twoWayRangeSearch(chrKey, start, end, Beds):
    values = []
    beds = {}
    for Bed in Beds:
        bed = getChr(Beds[Bed], chrKey)
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
            if intersect((chrKey, start, end), value)
        ]
    return values



def rangeMapping(toBed, fromBeds, DNAmet=False):
    chrs = set([bed[0] for bed in toBed.values()])
    if DNAmet:
        DNAmetMean = {}
        DNAmetLen = {}
    else:
        mappings = {}
        Proteins = list(fromBeds.keys())
    for chrKey in chrs:
        chrToBed = getChr(toBed, chrKey)
        chrToBed = sortBed(chrToBed)
        beds = dict.fromkeys(fromBeds.keys(), {})
        for fromBedKey in fromBeds:
            chrFromBed = getChr(fromBeds[fromBedKey], chrKey)
            chrFromBed = sortBed(chrFromBed)
            beds[fromBedKey] = chrFromBed
        for node in chrToBed:
            rangeX = chrToBed[node]
                nTuple = twoWayRangeSearch(rangeX[0], rangeX[1], rangeX[2], Beds=beds)
            if len(nTuple) != 0:
                # print(nTuple, rangeX)
                if DNAmet:
                    GCmetBed = fromBeds["met"]
                    C = [GCmetBed[met][3] for met, i in nTuple]
                    cov = [GCmetBed[met][4] for met, i in nTuple]
                    DNAmetMean[node] = (sum(C) / sum(cov))
                    DNAmetLen[node] = len(cov)
                else:
                    ProteinComb = [False for _ in Proteins]
                    for tup in nTuple:
                        idx = Proteins.index(tup[1])
                        ProteinComb[idx] = True
                    mappings[node] = ProteinComb
    if DNAmet:
        return DNAmetMean, DNAmetLen
    else:
        return Proteins, mappings





def oneWayRangeSearch(Bed, Protein, DNAmet=False):
    Proteins = tuple(Protein.keys())
    chrProteinStart = {}
    chrProteinEnd = {}
    chrProtein = {}
    values = {}
    idx = 0
    for protein, bed in Protein.items():
        for name, rangeX in bed.items():
            if rangeX[0] not in chrProtein:
                chrProteinStart[rangeX[0]] = []
                chrProteinEnd[rangeX[0]] = []
                chrProtein[rangeX[0]] = []
            chrProteinStart[rangeX[0]].append((rangeX[1], rangeX[2], protein, name))
            chrProteinEnd[rangeX[0]].append((rangeX[2], rangeX[1], protein, name))
            chrProtein[rangeX[0]].append((rangeX[0], rangeX[1], rangeX[2], protein, name))
    for key in chrProtein:
        chrProtein[key] = list(sorted(chrProtein[key], key=lambda v: v[1]))
        chrProteinEnd[key] = list(sorted(chrProteinEnd[key], key=lambda v: v[1]))
        chrProteinStart[key] = list(sorted(chrProteinStart[key], key=lambda v: v[0]))
    for cre, rangeX in Bed.items():
        idx += 1
        if idx % 1000 == 0:
            print(idx)
        left = bi.bisect_left(chrProteinEnd[rangeX[0]], (rangeX[1], rangeX[2]))
        right = bi.bisect_right(chrProteinStart[rangeX[0]], (rangeX[2], rangeX[1]))
        cre = tuple([cre] + list(rangeX))
        values[cre] = []
        values[cre] += [
            value
            for value in chrProtein[rangeX[0]][left:right]
            if intersect(rangeX, value[:3])
        ]
    mappings = {}
    for cre, bounds in values.items():
        boundComb = [False for _ in Proteins]
        creName = cre[0]
        for bound in bounds:
            boundName = bound[3]
            idx = Proteins.index(boundName)
            boundComb[idx] = True
        mappings[cre] = boundComb
    return values, mappings
