

def sortBed(Gbed):
    Gbed = {k: v for k,v in sorted(Gbed.items(), key=lambda kv: kv[1][1])}
    Gbed = {k: v for k, v in sorted(Gbed.items(), key=lambda kv: kv[1][0])}
    return Gbed


def writeBed(bed, file):
    f = open(file, "w+")
    for region in bed:
        f.write(bed[region][0]  + "\t" + str(bed[region][1]) + "\t" + str(bed[region][2]) + "\t" + region + "\n")
    f.close()


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

creMat = pd.read_csv("creCageConcat.tab", sep="\t")


tssCreP = {}
tssCreM = {}
mat = np.array(creMat.iloc[:,13:])
header = creMat.iloc[:,13:].columns
idxM = header.str.find("-") != -1
idxP = header.str.find("+") != -1
matM = mat[:,idxM]
matP = mat[:,idxP]



A = np.array([[1,2,3,4,5],
     [3,2,3,4,5],
     [4,2,3,4,5]])

A[:,::-1]

for creIdx, row in enumerate(matM):
    creChr = creMat.iloc[creIdx, 0]
    creStr = creMat.iloc[creIdx, 1]
    creEnd = creMat.iloc[creIdx, 2]
    creNam = creMat.iloc[creIdx, 3]
    tssCreM[(creChr, creStr, creEnd, creNam)] = []
    if creIdx % 1000 == 0:
        print(creIdx)
    for idx, cell in enumerate(row):
        head = header[idx]
        if row[idx-1] > cell and cell == 0:
            # idx is just after TSS, so idx-1 is the TSS
            # shift the center (500(left region) - idx(th bin at sample)*10(bin size))
            tssCreM[(creChr, creStr, creEnd, creNam)] += [(((creStr+creEnd)//2)-500) + (idx%100-1)*10]



matP = matP[:,::-1]

for creIdx, row in enumerate(matP):
    creChr = creMat.iloc[creIdx, 0]
    creStr = creMat.iloc[creIdx, 1]
    creEnd = creMat.iloc[creIdx, 2]
    creNam = creMat.iloc[creIdx, 3]
    tssCreP[(creChr, creStr, creEnd, creNam)] = []
    if creIdx % 1000 == 0:
        print(creIdx)
    for idx, cell in enumerate(row):
        head = header[idx]
        if row[idx-1] > cell and cell == 0:
            # idx is just after TSS, so idx-1 is the TSS
            # shift the center (500(left region) - idx(th bin at sample)*10(bin size))
            tssCreP[(creChr, creStr, creEnd, creNam)] += [(((creStr+creEnd)//2)+500) - (idx%100-1)*10]


ctrBedM = {}
ctrBedP = {}
ctrBed = {}
d = []
for key in tssCreP.keys():
    tssP = tssCreP[key]
    tssM = tssCreM[key]
    tmp = 99999999
    if len(tssP) == 0 and len(tssM) == 0:
        continue
    for i, p in enumerate(tssP):
        for j, m in enumerate(tssM):
            if tmp > abs(p - m):
                if abs(p - m) > 1000:
                    ex = ((i,p, tssP), (j,m, tssM), key)
                tmp = abs(p - m)
                tmpMin = (i, j)
    if len(tssP) == 0:
        tsM = max(tssCreM[key])
        ctrBedM[key[3]] = (key[0], tsM - 500, tsM + 500)
    elif len(tssM) == 0:
        tsP = min(tssCreP[key])
        ctrBedP[key[3]] = (key[0], tsP - 500, tsP + 500)
    else:
        i = tmpMin[0]
        j = tmpMin[1]
        tsP = (tssCreP[key][i])
        tsM = (tssCreM[key][j])
        # ctrBedP[key[3]] = (key[0], tsP - 500, tsP + 500)
        # ctrBedM[key[3]] = (key[0], tsM - 500, tsM + 500)
        if tsP > tsM:
            ctrBed[key[3]] = (key[0], tsM, tsP)
        elif tsP == tsM:
            ctrBed[key[3]] = (key[0], tsP, tsM+2)
        else:
            ctrBed[key[3]] = (key[0], tsP, tsM)
        d += [tmp]




ctrBedM = sortBed(ctrBedM)
writeBed(ctrBedM, "creCtrM.bed")

ctrBedP = sortBed(ctrBedP)
writeBed(ctrBedP, "creCtrP.bed")

ctrBed = sortBed(ctrBed)
writeBed(ctrBed, "creCtr.bed")





Z = d
N = len(Z)
X2 = np.sort(Z)
F2 = np.array(range(N))/float(N)
plt.plot(X2, F2, label="cre" + ": " + str(N), color=".25")


plt.legend()

plt.ylabel("CDF")
plt.xlabel("D")
plt.title("Distance")
plt.show()
