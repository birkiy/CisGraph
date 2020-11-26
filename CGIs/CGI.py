

from Functions.Packages import *




conAnn = pd.read_csv(f"{dataRoot}/Regions/con.AnnFile.cvs", sep="\t")
indAnn = pd.read_csv(f"{dataRoot}/Regions/ind.AnnFile.csv", sep="\t")
nonAnn = pd.read_csv(f"{dataRoot}/Regions/non.AnnFile.csv", sep="\t")
othAnn = pd.read_csv(f"{dataRoot}/Regions/oth.AnnFile.csv", sep="\t")
tsPAnn = pd.read_csv("/home/ualtintas/genomeAnnotations/Regions/TSS.hg19.+.AnnFile.csv", sep="\t")
tsMAnn = pd.read_csv("/home/ualtintas/genomeAnnotations/Regions/TSS.hg19.-.AnnFile.csv", sep="\t")

conAnn = conAnn.rename(columns={"seqnames": "Chr", "start": "Start", "end": "End"})
indAnn = indAnn.rename(columns={"seqnames": "Chr", "start": "Start", "end": "End"})
nonAnn = nonAnn.rename(columns={"seqnames": "Chr", "start": "Start", "end": "End"})
othAnn = othAnn.rename(columns={"seqnames": "Chr", "start": "Start", "end": "End"})
tsPAnn = tsPAnn.rename(columns={"seqnames": "Chr", "start": "Start", "end": "End"})
tsMAnn = tsMAnn.rename(columns={"seqnames": "Chr", "start": "Start", "end": "End"})



con = conAnn[["V4", "Chr", "Start", "End", "annotation", "geneStrand", "SYMBOL"]]
ind = indAnn[["V4", "Chr", "Start", "End", "annotation", "geneStrand", "SYMBOL"]]
non = nonAnn[["V4", "Chr", "Start", "End", "annotation", "geneStrand", "SYMBOL"]]
oth = othAnn[["V4", "Chr", "Start", "End", "annotation", "geneStrand", "SYMBOL"]]
tsP = tsPAnn[["V4", "Chr", "Start", "End", "annotation", "geneStrand", "SYMBOL"]]
tsM = tsMAnn[["V4", "Chr", "Start", "End", "annotation", "geneStrand", "SYMBOL"]]


con["nodeClass"] = "con"
ind["nodeClass"] = "ind"
non["nodeClass"] = "non"
oth["nodeClass"] = "oth"
tsP["nodeClass"] = "tsP"
tsM["nodeClass"] = "tsM"


# allAnn = pd.concat([con, ind, non, oth, tsP, tsM]). reset_index(drop=True)
allAnn = pd.concat([con, ind, non]). reset_index(drop=True)


allAnn["annotation"] = allAnn["annotation"].str.split("(", expand=True)[0]
allAnn = allAnn.sort_values("Start").reset_index(drop=True)


cpg = pd.read_table("/home/ualtintas/genomeAnnotations/cpgIslandExt.bed", sep="\t", names=["Chr", "Start", "End", "Name", "Length", "cpgNum", "gcNum", "perCpg", "perGc", "obsExp"])
cpg = cpg.sort_values("Start").reset_index(drop=True)





dictCpG = {}
for i in range(len(allAnn)):
    chr = allAnn.loc[i, "Chr"]
    start = allAnn.loc[i, "Start"]
    end = allAnn.loc[i, "End"]
    leftV = bi.bisect_right(cpg["Start"] , end)
    rightV = bi.bisect_left(cpg["End"] , start)
    tmp = cpg[rightV:leftV]
    tmp = tmp[tmp["Chr"] == chr]
    nodeClass = allAnn.loc[i, "nodeClass"]
    if nodeClass not in dictCpG:
        dictCpG[nodeClass] = tmp.shape[0]
        if tmp.shape[0] != 0:
            print(tmp)
            print(chr, start, end)
    else:
        dictCpG[nodeClass] += tmp.shape[0]
        if tmp.shape[0] != 0:
            print(tmp)
            print(chr, start, end)
    if i % 1000 == 0:
        print(dictCpG)



    
