
from Functions.Packages import *


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

def getChr(Gbed, chrKey):
    Gbed = {k: v for k, v in Gbed.items() if v[0] == chrKey}
    return Gbed

def nodeMapping(toBed, fromBeds):
    chrs = set([bed[0] for bed in toBed.values()])
    mappings = {}
    for chrKey in chrs:
        chrToBed = getChr(toBed, chrKey)
        beds = dict.fromkeys(formBeds.keys(), {})
        for fromBedKey in fromBeds:
            chrFromBed = getChr(fromBeds[fromBedKey], chrKey)
            beds[formBedKey] = chrFromBed
        for node in chrToBed:
            rangeX = chrToBed[node]
            nTuple = rangesFromUpperRange(rangeX[0], rangeX[1], rangeX[2], Beds=beds)
            if len(nTuple) != 0:
                print(nTuple, rangeX)
                mappings[node] = nTuple
    return mappings



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
