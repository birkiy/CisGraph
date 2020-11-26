


from Functions.Packages import *


rip = pd.read_csv(f"{dataRoot}/GSE100710_RIP-Seq_DeSEQ2_AR_x_IgG.tab.txt", sep="\t")


ripDE = rip[((rip["RIP_log2(AR/IgG)"] > 1) | (rip["RIP_log2(AR/IgG)"] < -1)) & (rip["RIP_log2(AR/IgG)_adjusted_p_value"] < 0.05)]


new = ripDE["locus"].str.split(":", expand=True)
ripDE["chr1"] = new[0]

new1 = new[1].str.split("-", expand=True)
ripDE["start"] = new1[0]
ripDE["end"] = new1[1]

ripDE["chr2"] = "chr"
ripDE["chr"] = ripDE["chr2"] + ripDE["chr1"].astype(str)
ripDE = ripDE.drop(["chr1", "chr2"], axis=1)

ripDE = ripDE.astype({"start": "int32", "end": "int32"})
ripDE = ripDE.sort_values("start").reset_index(drop=True)



conBed = pd.read_csv(home + "/ARBSs/regions/cons-arbs.bed", sep="\t", names=["chr", "start", "end", "name"])
conBed["nodeClass"] = "con"
indBed = pd.read_csv(home + "/ARBSs/regions/ind-arbs.bed", sep="\t", names=["chr", "start", "end", "name"])
indBed["nodeClass"] = "ind"
nonBed = pd.read_csv(home + "/ARBSs/regions/Non-Active-ARBS.bed", sep="\t", names=["chr", "start", "end", "name"])
nonBed["nodeClass"] = "non"
nARBed = pd.read_csv(home + "/ARBSs/regions/negativeControl.ARBS.bed", sep="\t", names=["chr", "start", "end", "name"])
nARBed["nodeClass"] = "nAR"

tsPBed = pd.read_csv(home + "/genomeAnnotations/Regions/TSS.hg19.+.bed", sep="\t", names=["chr", "start", "end", "name", "strand"])
tsPBed["nodeClass"] = "tss"
tsMBed = pd.read_csv(home + "/genomeAnnotations/Regions/TSS.hg19.-.bed", sep="\t", names=["chr", "start", "end", "name", "strand"])
tsMBed["nodeClass"] = "tss"


DF = pd.concat([conBed, indBed, nonBed, nARBed, tsPBed, tsMBed])
DF["center"] = (DF["start"] + DF["end"]) // 2

DF = DF.sort_values("center").reset_index(drop=True)

pickle.dump(DF, open(f"{dataRoot}/tmpData/BED.DF.p", "wb"))


sel = {"-": "end", "+": "start"}

map = {}
for i in range(len(ripDE)):
    chrKey = ripDE.loc[i, "chr"]
    if chrKey not in map:
        map[chrKey] = []
    s = ripDE.loc[i,"strand"]
    col = sel[s]
    if s == "-":
        map[chrKey].append(ripDE.loc[i, "end"])
    elif s == "+":
        map[chrKey].append(ripDE.loc[i, "start"])


def find_le(a, x):
    'Find rightmost value less than or equal to x'
    i = bi.bisect_right(a, x)
    if i:
        return a[i-1]
    raise ValueError

def find_ge(a, x):
    'Find leftmost item greater than or equal to x'
    i = bi.bisect_left(a, x)
    if i != len(a):
        return a[i]
    raise ValueError


DF["minDis"] = 999999999
for i in range(len(DF)):
    center = DF.loc[i, "center"]
    cChr = DF.loc[i, "chr"]
    cStrand = DF.loc[i, "strand"]
    if cChr not in map:
        continue
    minVal = 10000
    try:
        right = find_ge(map[cChr], center)
        minVal = min(minVal, abs(center - right))
    except:
        pass
    try:
        left = find_le(map[cChr], center)
        minVal = min(minVal, abs(center - left))
    except:
        pass
    DF.loc[i,"minDis"] = minVal
    if i%5000 == 0:
        print(i)




DF = DF[DF["minDis"] < 10000]


fig = plt.figure(figsize=[5,6])

gs = gridspec.GridSpec(ncols=1, nrows=1)

colorPalette = {"con": "#63b7af",
                "ind": "#abf0e9",
                "non": "#d4f3ef",
                "nAR": "#f5fffd",
                "tss": "#ee8572"}


fig.add_subplot(gs[0])
for nodeClass, color in colorPalette.items():
    Z = DF[DF["nodeClass"] == nodeClass]["minDis"]
    N = len(Z)
    X2 = np.sort(Z)
    F2 = np.array(range(N))/float(N)
    plt.plot(X2, F2, label=nodeClass + ": " + str(N), color=color)


plt.legend()

plt.ylabel("CDF")
plt.xlabel("Distance")

fig.savefig(f"{figureRoot}/CDF.pdf")










q = [1,3,5]
r = bi.bisect_left(q, 4)
q[l]

l = bi.bisect_right(q, 4) -1
q[r]
