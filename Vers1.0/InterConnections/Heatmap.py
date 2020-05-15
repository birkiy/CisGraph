


from Functions.PlotFunctions import *
from Functions.Packages import *

home = "/home/birkiy/github/CisGraph"

R = pickle.load(open(home + "/Data/tmpData/InterT.RData.p", "rb" ))
S = pickle.load(open(home + "/Data/tmpData/InterT.SData.p", "rb" ))


fig = plt.figure()
fig.subplots_adjust(hspace=0.2, wspace=0.2)
ax = fig.add_subplot(1, 1, 1)

hc = ["#2c6589", "#87a7b9", "#ffffff", "#aa5f55", "#95342c"]
th = [0, 0.15, 0.5, 0.85, 1]
cdict = NonLinCdict(th, hc)
cm = colors.LinearSegmentedColormap('test', cdict)

sns.heatmap(R, xticklabels=["hom", "con", "ind", "non"],
            yticklabels=["hom", "con", "ind", "non"],
            cmap=cm, center=1)

fig.savefig(home + "/Figures/InterTADHeatmap.pdf")

plt.close("all")




fig = plt.figure()
fig.subplots_adjust(hspace=0.2, wspace=0.2)
ax = fig.add_subplot(1, 1, 1)


sns.heatmap(S, xticklabels=["hom", "con", "ind", "non"],
            yticklabels=["hom", "con", "ind", "non"],
            cmap="Greys", center=1)

fig.savefig(home + "/Figures/InterTADHeatmapS.pdf")

plt.close("all")
