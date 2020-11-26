

from Functions.Packages import *

chipARBS = pd.read_csv(f"{dataRoot}/ARBS_chip_seq_2.csv")



biChipARBS = chipAll.set_index("Unnamed: 0")



for col in biChipARBS:
    biChipARBS[col] = biChipARBS.apply(lambda b: 1 if b[col] > 0.5 else 0, axis=1)


biChipARBS.sum(axis=1)

nodeClass = mat[:,-1]

mat = biChipARBS.values
motifAR = {}
motifNode = {"con": {}, "ind": {}, "non": {}, "nAR": {}, "tss": {}}
for i, row in enumerate(mat):
    row = tuple(row[:-1])
    if row not in motifAR:
        motifAR[row] = []
    if row not in motifNode[mat[i][-1]]:
        motifNode[mat[i][-1]][row] = []
    motifAR[row].append(i)
    motifNode[mat[i][-1]][row].append(i)


[len(x) for x in motifNode["con"].values()]

check = biChipARBS.shape[1]

motifAR = {}
visited = []
c = 0
remain = list(range(biChipARBS.shape[0]))
for row in range(biChipARBS.shape[0]):
    motifAR[row] = [row]
    if row in visited:
        continue
    visited += [row]
    remain.remove(row)
    tmp = remain
    if remain == []:
        break
    for row2 in remain:
        if sum(mat[row] == mat[row2]) == check:
            motifAR[row] += [row2]
            visited += [row2]
            tmp.remove(row2)
    remain = tmp
    if c % 1000 == 0:
        print(c)
    c += 1


cpy = motifAR.copy()
for key, value in cpy.items():
    if value == []:
        del motifAR[key]





len(motifAR.keys()) + sum([len(x) for x in motifAR.values()])


biChipARBS["type"] = -1
for key, value in motifAR.items():
    idx = biChipARBS.index[key]
    biChipARBS.loc[idx, "type"] = key
    for i in value:
        idx = biChipARBS.index[i]
        biChipARBS.loc[idx, "type"] = key


DF = biChipARBS.groupby(["nodeClass", "type"]).size()

DF = DF.reset_index()
DF.sort_values(0)

DF[(DF[0] != 1) & (DF["type"] != -1 )].sort_values(0)

pd.set_option("display.max_rows", None)
