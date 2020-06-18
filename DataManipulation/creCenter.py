

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
            tssCreM[(creChr, creStr, creEnd, creNam)] += [((((creStr+creEnd)//2)-500) + (idx%100-1)*10, row[idx-1])]



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
            tssCreP[(creChr, creStr, creEnd, creNam)] += [((((creStr+creEnd)//2)+500) - (idx%100-1)*10, row[idx-1])]


ctrBedM = {}
ctrBedM1 = {}
ctrBedM2 = {}

ctrBedP = {}
ctrBedP1 = {}
ctrBedP2 = {}

ctrBedA = {}
ctrBedA1 = {}
ctrBedA2 = {}

ctrBedB = {}
ctrBedB1 = {}
ctrBedB2 = {}

d = []
for key in tssCreP.keys():
    tssP = tssCreP[key]
    tssM = tssCreM[key]
    tmpPE = 0
    tmpME = 0
    for i, p in enumerate(tssP):
        pP, pE = p
        if pE >= tmpPE:
            tmpPE = pE
            tsP = pP
    for j, m in enumerate(tssM):
        mP, mE = m
        if mE >= tmpME:
            tmpME = mE
            tsM = mP
    if len(tssP) == 0 and len(tssM) == 0:
        continue
    if len(tssP) == 0:
        ctrBedM[key[3]] = (key[0], tsM - 500, tsM + 500)
        ctrBedM1[key[3]] = (key[0], tsM - 500, tsM)
        ctrBedM2[key[3]] = (key[0], tsM, tsM + 500)
    elif len(tssM) == 0:
        ctrBedP[key[3]] = (key[0], tsP - 500, tsP + 500)
        ctrBedP1[key[3]] = (key[0], tsP - 500, tsP)
        ctrBedP2[key[3]] = (key[0], tsP, tsP + 500)
    else:
        tmp = abs(tmpME - tmpPE)
        # ctrBedP[key[3]] = (key[0], tsP - 500, tsP + 500)
        # ctrBedM[key[3]] = (key[0], tsM - 500, tsM + 500)
        if tsP >= tsM:
            ctrBedA[key[3]] = (key[0], tsM-400, tsP+400)
            ctrBedA1[key[3]] = (key[0], tsM-400, (tsP+tsM)//2)
            ctrBedA2[key[3]] = (key[0], (tsP+tsM)//2, tsP+400)
        else:
            ctrBedB[key[3]] = (key[0], tsP-400, tsM+400)
            ctrBedB1[key[3]] = (key[0], tsP-400, (tsP+tsM)//2)
            ctrBedB2[key[3]] = (key[0], (tsP+tsM)//2, tsM+400)
        d += [tmp]




ctrBedM = sortBed(ctrBedM)
writeBed(ctrBedM, "creCtrM.bed")

ctrBedM1 = sortBed(ctrBedM1)
writeBed(ctrBedM1, "creCtrM1.bed")

ctrBedM2 = sortBed(ctrBedM2)
writeBed(ctrBedM2, "creCtrM2.bed")



ctrBedP = sortBed(ctrBedP)
writeBed(ctrBedP, "creCtrP.bed")

ctrBedP1 = sortBed(ctrBedP1)
writeBed(ctrBedP1, "creCtrP1.bed")

ctrBedP2 = sortBed(ctrBedP2)
writeBed(ctrBedP2, "creCtrP2.bed")




ctrBedA = sortBed(ctrBedA)
writeBed(ctrBedA, "creCtrA.bed")

ctrBedA1 = sortBed(ctrBedA1)
writeBed(ctrBedA1, "creCtrA1.bed")

ctrBedA2 = sortBed(ctrBedA2)
writeBed(ctrBedA2, "creCtrA2.bed")



ctrBedB = sortBed(ctrBedB)
writeBed(ctrBedB, "creCtrB.bed")

ctrBedB1 = sortBed(ctrBedB1)
writeBed(ctrBedB1, "creCtrB1.bed")

ctrBedB2 = sortBed(ctrBedB2)
writeBed(ctrBedB2, "creCtrB2.bed")

creMat.iloc[:,[3,12]]

ctrBed = {**ctrBedB, **ctrBedA, **ctrBedM, **ctrBedP}
ctrBed = sortBed(ctrBed)
writeBed(ctrBed, "creCtr.bed")


creDHT = {k: v for k, v in ctrBed.items() if k.find("DHT") != -1}
writeBed(creDHT, "creDHT.bed")
creETH = {k: v for k, v in ctrBed.items() if k.find("EtOH") != -1}
writeBed(creETH, "creEtOH.bed")


Z = d
N = len(Z)
X2 = np.sort(Z)
F2 = np.array(range(N))/float(N)
plt.plot(X2, F2, label="cre" + ": " + str(N), color=".25")


plt.legend()

plt.ylabel("CDF")
plt.xlabel("D")
plt.title("Distance")

plt.savefig("distance.pdf")
