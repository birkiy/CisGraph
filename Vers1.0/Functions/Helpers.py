
from Functions.Packages import *


def readBed(bed, file):
    with open(file) as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        for row in reader:
            # row3 = row[3].split(".")[0]
            row3 = row[3]
            bed[row3] = (row[0], int(row[1]), int(row[2]))



def rangesFromUpperRange(chrKey, start, end, beds, bedname):
    values = []
    for idx, enhancer in enumerate(beds):
        enhancerEnd = [(value[2], value[1], value[0], key) for key, value in enhancer.items()]
        left = bi.bisect_left(enhancerEnd, (start, end))

        enhancerStart = [(value[1], value[2], value[0], key) for key, value in enhancer.items()]
        right = bi.bisect_right(enhancerStart, (end, start))
        values += [
            (keys, bedname[idx])
            for keys, value in list(enhancer.items())[left:right]
            if chrKey == value[0] and
               value[2] - start > start - value[1] and
               end - value[1] > value[2] - end
        ]
    return values
