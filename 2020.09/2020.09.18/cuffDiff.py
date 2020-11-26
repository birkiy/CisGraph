

from Functions.Packages import *

ripAR = pd.read_csv("/home/ualtintas/STARexample/diff.out/isoform_exp.diff", sep="\t")

"""
CuffDiff isoforms output gives the transcripts not the genes, meaning alternative splicing,
but still I could not get the same results with suplementary data. I continued with supplementary data.
"""


ripDE = ripAR[(ripAR["sample_1"] == "AR") & (ripAR["sample_2"] == "IgG") & (ripAR["status"] == "OK") & (ripAR["significant"] == "yes")].sort_values("value_1").reset_index(drop=True)


ripAR[(ripAR["sample_1"] == "AR") & (ripAR["sample_2"] == "IgG") & (ripAR["status"] == "OK") & (ripAR["q_value"] < 0.05) & ((ripAR["log2(fold_change)"] > 1) | (ripAR["log2(fold_change)"] < -1))]


new = ripDE["locus"].str.split(":", expand=True)

ripDE["chr"] = new[0]

new2 = new[1].str.split("-", expand=True)
ripDE["start"] = new2[0]
ripDE["end"] = new2[1]

ripDE = ripDE.sort_values("start").reset_index(drop=True)

bed = pickle.load(open(f"{dataRoot}/tmpData/BED.DF.p", "rb"))




map = {}
for i in range(len(ripDE)):
    chrKey = ripDE.loc[i, "chr"]
    if chrKey not in map:
        map[chrKey] = []
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


bed["minDis"] = 999999999
for i in range(len(bed)):
    center = bed.loc[i, "center"]
    cChr = bed.loc[i, "chr"]
    cStrand = bed.loc[i, "strand"]
    if cChr not in map:
        continue
    minVal = 1000000
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
    bed.loc[i,"minDis"] = minVal
    if i%5000 == 0:
        print(i)

minBed = bed[bed["minDis"] < 1000000]


fig = plt.figure(figsize=[4, 5])
gs = gridspec.GridSpec(ncols=1, nrows=1)
plt.subplots_adjust(wspace=0.4)

colorPalette = ["#63b7af", "#abf0e9", "#d4f3ef", "#f5fffd", "#ee8572"]

params = dict(x='nodeClass',
              y='minDis',
              order=["con", "ind", "non", "nAR", "tss"])

boxPairs = [("nAR", "con"),
            ("con", "ind"),
            ("con", "non"),
            ("ind", "non"),
            ("con", "tss"),
            ("non", "tss"),
            ("nAR", "tss")
            ]

ax1 = fig.add_subplot(gs[0])
ax1 = sns.boxplot(**params, data=minBed, palette=colorPalette)
plt.title("Minimum distance to \n AR bound transcript code")
plt.ylabel("Min distance (bp)", fontsize=14)
plt.xlabel("Node Classes", fontsize=14)
add_stat_annotation(ax1, **params, data=minBed, box_pairs=boxPairs, test='Mann-Whitney')


fig.savefig(f"{figureRoot}/minDis.pdf")
