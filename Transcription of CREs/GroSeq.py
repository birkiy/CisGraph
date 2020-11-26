

from Functions.Packages import *

con = pd.read_csv(f"{dataRoot}/Gro/cons-arbs.tab", sep="\t", names=["dmso.-","dmso.+","dht.-","dht.+"])
con["nodeClass"] = "con"
ind = pd.read_csv(f"{dataRoot}/Gro/ind-arbs.tab", sep="\t", names=["dmso.-","dmso.+","dht.-","dht.+"])
ind["nodeClass"] = "ind"
nAR = pd.read_csv(f"{dataRoot}/Gro/negativeControl.ARBS.tab", sep="\t", names=["dmso.-","dmso.+","dht.-","dht.+"])
nAR["nodeClass"] = "nAR"
non = pd.read_csv(f"{dataRoot}/Gro/Non-Active-ARBS.tab", sep="\t", names=["dmso.-","dmso.+","dht.-","dht.+"])
non["nodeClass"] = "non"
tss = pd.read_csv(f"{dataRoot}/Gro/TSS.hg19.Idx.woS.tab", sep="\t", names=["dmso.-","dmso.+","dht.-","dht.+"])
tss["nodeClass"] = "pro"


DF = pd.concat([con, ind, non, nAR, tss])

pickle.dump(DF, open(f"{dataRoot}/tmpData/Gro.DF.p", "wb"))


GroDF = pickle.load(open(f"{dataRoot}/tmpData/Gro.DF.p", "rb" ))

GroDF = GroDF.reset_index()
GroDF = GroDF.rename(columns={"index": "Gene"})
GroDF = GroDF.sort_values("Gene")
GroDF = GroDF.reset_index().drop("index", axis=1)

GroDF["dmso"] = np.log(GroDF["dmso.-"] + GroDF["dmso.+"] + 1)

GroDF["dht"] = np.log(GroDF["dht.-"] + GroDF["dht.+"] + 1)



#############################################################

fig = plt.figure(figsize=[8, 5])
gs = gridspec.GridSpec(ncols=2, nrows=1)
plt.subplots_adjust(wspace=0.4)


#colorPalette = ["#63b7af", "#abf0e9","#d4f3ef", "#f5fffd", "#ee8572"]

colorPalette = {"con": "#5A5A5A",
                "ind": "#F9746D",
                "non": "#ACACAC",
                "nAR": "#F5F5F5",
                "pro": "#000000"}

params = dict(x='nodeClass',
              order=["con", "ind", "non", "nAR", "pro"],
              data=GroDF)


boxPairs = [("nAR", "con"),
            ("con", "ind"),
            ("con", "non"),
            ("ind", "non"),
            ("con", "pro"),
            ("non", "pro"),
            ("nAR", "pro")
            ]


ax1 = fig.add_subplot(gs[0])
ax1 = sns.boxplot(**params, y="dmso", palette=colorPalette)
plt.title("DMSO")
plt.ylabel("GroSeq Signal Level (Log Scaled)", fontsize=14)
plt.xlabel("Node Classes", fontsize=14)
for tick in ax1.xaxis.get_major_ticks():
    tick.label.set_fontsize(13)


add_stat_annotation(ax1, **params, y="dmso", box_pairs=boxPairs, test='Mann-Whitney')

ax2 = fig.add_subplot(gs[1])
ax2 = sns.boxplot(**params, y="dht", palette=colorPalette)
plt.title("DHT")
plt.ylabel("GroSeq Signal Level (Log Scaled)", fontsize=14)
plt.xlabel("Node Classes", fontsize=14)
for tick in ax2.xaxis.get_major_ticks():
    tick.label.set_fontsize(13)


add_stat_annotation(ax2, **params, y="dmso", box_pairs=boxPairs, test='Mann-Whitney')

fig.savefig(f"{figureRoot}/Gro.pdf")
