

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

creMat = pd.read_csv("creCageConcat.tab", sep="\t")

#
# samples = ["LNCaP.0h.rep1", "LNCaP.3h.rep1", "LNCaP.6h.rep1", "LNCaP.12h.rep1", "LNCaP.18h.rep1", "LNCaP.48h.rep1", "LNCaP.72h.rep1"]
# extends = [".-.BPM", ".+.BPM"]
#
# for sample in samples:
#     for extend in extends:
#         for i in range(100):
#             col = sample + extend

tssCreP = {}
tssCreM = {}
mat = np.array(creMat.iloc[:,13:])
header = creMat.iloc[:,13:].columns
for creIdx, row in enumerate(mat):
    creChr = creMat.iloc[creIdx, 0]
    creStr = creMat.iloc[creIdx, 1]
    creEnd = creMat.iloc[creIdx, 2]
    creNam = creMat.iloc[creIdx, 3]
    tssCreP[(creChr, creStr, creEnd, creNam)] = []
    tssCreM[(creChr, creStr, creEnd, creNam)] = []
    if creIdx % 1000 == 0:
        print(creIdx)
    for idx, cell in enumerate(row):
        head = header[idx]
        # if idx % 100 == 0:# print("\t" + str(idx))
        if head.find("-") == -1:
            if row[idx-1] < cell and row[idx-1] == 0:
                # idx is TSS
                tssCreM[(creChr, creStr, creEnd, creNam)] += [(((creStr+creEnd)//2)-500) + idx*10]
        if head.find("+") == -1:
            if row[idx-1] > cell and cell == 0:
                #idx is TSS
                tssCreP[(creChr, creStr, creEnd, creNam)] += [(((creStr+creEnd)//2)-500) + (idx-1)*10]



# tssCrePsorted = dict(sorted(tssCreP.items(), key= lambda kv: kv[0][3]))
# tssCreMsorted = dict(sorted(tssCreM.items(), key= lambda kv: kv[0][3]))

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
                tmp = abs(p - m)
                tmpMin = (i, j)
    if len(tssP) == 0:
        tsM = min(tssCreM[key])
        ctrBedM[key[3]] = (key[0], tsM - 500, tsM + 500)
    elif len(tssM) == 0:
        tsP = max(tssCreP[key])
        ctrBedP[key[3]] = (key[0], tsP - 500, tsP + 500)
    else:
        i = tmpMin[0]
        j = tmpMin[1]
        tsP = (tssCreP[key][i])
        tsM = (tssCreM[key][j])
        ctrBedP[key[3]] = (key[0], tsP - 500, tsP + 500)
        ctrBedM[key[3]] = (key[0], tsM - 500, tsM + 500)
        if tsP > tsM:
            ctrBed[key[3]] = (key[0], tsM, tsP)
        elif tsP == tsM:
            ctrBed[key[3]] = (key[0], tsP, tsM+1)
        else:
            ctrBed[key[3]] = (key[0], tsP, tsM)
        d += [tmp]




ctrBedM = sortBed(ctrBedM)
writeBed(ctrBedM, "creCtrM.bed")

ctrBedP = sortBed(ctrBedP)
writeBed(ctrBedP, "creCtrP.bed")

ctrBed = sortBed(ctrBed)
writeBed(ctrBed, "creCtr.bed")




## MAX
ctrBedM = {}
ctrBedP = {}
ctrBed = {}
d = []
for key in tssCreP.keys():
    tssP = tssCreP[key]
    tssM = tssCreM[key]
    tmp = 0
    if len(tssP) == 0 and len(tssM) == 0:
        continue
    for i, p in enumerate(tssP):
        for j, m in enumerate(tssM):
            if tmp < abs(p - m):
                tmp = abs(p - m)
                tmpMax = (i, j)
    if len(tssP) == 0:
        tsM = min(tssCreM[key])
        ctrBedM[key[3]] = (key[0], tsM - 500, tsM + 500)
    elif len(tssM) == 0:
        tsP = max(tssCreP[key])
        ctrBedP[key[3]] = (key[0], tsP - 500, tsP + 500)
    else:
        i = tmpMax[0]
        j = tmpMax[1]
        tsP = (tssCreP[key][i])
        tsM = (tssCreM[key][j])
        ctrBedP[key[3]] = (key[0], tsP - 500, tsP + 500)
        ctrBedM[key[3]] = (key[0], tsM - 500, tsM + 500)
        if tsP > tsM:
            ctrBed[key[3]] = (key[0], tsM, tsP)
        elif tsP == tsM:
            ctrBed[key[3]] = (key[0], tsP, tsM+1)
        else:
            ctrBed[key[3]] = (key[0], tsP, tsM)
        d += [tmp]

ctrBedM = sortBed(ctrBedM)
writeBed(ctrBedM, "creCtrM.bed")

ctrBedP = sortBed(ctrBedP)
writeBed(ctrBedP, "creCtrP.bed")

ctrBed = sortBed(ctrBed)
writeBed(ctrBed, "creCtr.bed")



# x=np.random.normal(10, 0.1, size=(5000))

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
