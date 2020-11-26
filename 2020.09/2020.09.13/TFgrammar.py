
from Functions.Packages import *


header = ["chr", "start", "end", "NumOfBound", "Which", "AR", "ARID1A", "CHD1", "CHD4", "CSNK2A1", "CTBP1", "CTBP2", "CTCF", "EZH2", "FOXA1", "GATA2", "GRHL2", "HIF1A", "HOXB13", "KDM1A", "MED", "MTOR", "NKX31", "NKX3", "NMYC", "PIAS1", "POU2F1", "RELA", "RNAP2", "SFPQ", "SMARCA1", "SMARCA2", "SMARCA4", "SMARCA5", "SUZ12", "TCF7L2", "TET2", "TFAP4", "TLE3", "TRIM24", "TRIM28", "WDHD1", "WDR5", "ARBS.chr", "ARBS.start", "ARBS.end", "ARBS.name"]
TFs = ["AR", "ARID1A", "CHD1", "CHD4", "CSNK2A1", "CTBP1", "CTBP2", "CTCF", "EZH2", "FOXA1", "GATA2", "GRHL2", "HIF1A", "HOXB13", "KDM1A", "MED", "MTOR", "NKX31", "NKX3", "NMYC", "PIAS1", "POU2F1", "RELA", "RNAP2", "SFPQ", "SMARCA1", "SMARCA2", "SMARCA4", "SMARCA5", "SUZ12", "TCF7L2", "TET2", "TFAP4", "TLE3", "TRIM24", "TRIM28", "WDHD1", "WDR5"]


con = pd.read_csv(home + "/ARBSs/cre/con.cre.bed", sep="\t", names=header)
ind = pd.read_csv(home + "/ARBSs/cre/ind.cre.bed", sep="\t", names=header)
nAR = pd.read_csv(home + "/ARBSs/cre/nAR.cre.bed", sep="\t", names=header)
non = pd.read_csv(home + "/ARBSs/cre/non.cre.bed", sep="\t", names=header)
tss = pd.read_csv(home + "/ARBSs/cre/tss.cre.bed", sep="\t", names=header + ["strand"])




biCon = con.groupby("ARBS.name").sum()
biInd = ind.groupby("ARBS.name").sum()
biNon = non.groupby("ARBS.name").sum()
biNAR = nAR.groupby("ARBS.name").sum()
biTss = tss.groupby("ARBS.name").sum()


biCon = biCon.drop(columns=["start", "end", "ARBS.start", "ARBS.end"])
biInd = biInd.drop(columns=["start", "end", "ARBS.start", "ARBS.end"])
biNon = biNon.drop(columns=["start", "end", "ARBS.start", "ARBS.end"])
biNAR = biNAR.drop(columns=["start", "end", "ARBS.start", "ARBS.end"])
biTss = biTss.drop(columns=["start", "end", "ARBS.start", "ARBS.end"])


for col in TFs:
    biCon[col] = biCon.apply(lambda b: 0 if b[col] == 0 else 1, axis=1)
    biInd[col] = biInd.apply(lambda b: 0 if b[col] == 0 else 1, axis=1)
    biNon[col] = biNon.apply(lambda b: 0 if b[col] == 0 else 1, axis=1)
    biNAR[col] = biNAR.apply(lambda b: 0 if b[col] == 0 else 1, axis=1)
    biTss[col] = biTss.apply(lambda b: 0 if b[col] == 0 else 1, axis=1)



biCon["nodeClass"] = "con"
biInd["nodeClass"] = "ind"
biNon["nodeClass"] = "non"
biNAR["nodeClass"] = "nAR"
biTss["nodeClass"] = "tss"



biCRE = pd.concat([biCon, biInd, biNon, biNAR, biTss])

##########################################

nodeClass = mat[:,-1]

mat = biCRE.values
motifAR = {}
motifNode = {}
for i, row in enumerate(mat):
    row = tuple(row[1:-1])
    if row not in motifAR:
        motifAR[row] = []
    if row not in motifNode:
        motifNode[row] = {"con": [], "ind": [], "non": [], "nAR": [], "tss": []}
    motifAR[row].append(i)
    motifNode[row][mat[i][-1]].append(i)


cpyMotif = motifNode.copy()
for motif in cpyMotif.keys():
    if motifNode[motif]["con"] == [] and motifNode[motif]["ind"] == [] and motifNode[motif]["non"] == []:
        del motifNode[motif]





cats = [[len(x) for x in motifDict.values()] for motifDict in motifNode.values()]

cpyMotif = motifNode.copy()
for i, motif in enumerate(cpyMotif.keys()):
    if sum(cats[i]) == 1:
        del motifNode[motif]



TFs = ["AR", "PIAS1", "ARID1A", "MED", "SMARCA4", "CTBP2", "TLE3", "RNAP2", "GATA2", "WDHD1", "FOXA1"]
TFs = ["AR"]

smCRE = biCRE[TFs + ["nodeClass"]]
mat = smCRE.values

motifNode2 = {}
for i, row in enumerate(mat):
    row = tuple(row[:-1])
    if row not in motifNode2:
        motifNode2[row] = {"con": [], "ind": [], "non": [], "nAR": [], "tss": []}
    motifNode2[row][mat[i][-1]].append(i)

[[len(x) for x in motifDict.values()] for motifDict in motifNode2.values()]


MRE11A
RUNX1
