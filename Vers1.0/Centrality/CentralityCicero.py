
from Functions.Packages import *




ethG = pickle.load(open(f"{dataRoot}/tmpData/GraphsG.EtOH.Data.p", "rb" ))
dhtG = pickle.load(open(f"{dataRoot}/tmpData/GraphsG.DHT.Data.p", "rb" ))




# Not Connected
ethb = nx.betweenness_centrality(ethG)
ethd = nx.degree_centrality(ethG)

ethB = pd.DataFrame.from_dict(ethb, orient="index", columns=["Betweenness Centrality"])

ethD = pd.DataFrame.from_dict(ethd, orient="index", columns=["Degree Centrality"])

ethB["Degree Centrality"] = ethD["Degree Centrality"]

nodeClasses = []
for node in ethB.index:
    nodeClasses += [ethG.nodes()[node]["nodeClass"]]


ethB["nodeClass"] = nodeClasses
ethB["Condition"] = "EtOH"


# Not Connected
dhtb = nx.betweenness_centrality(dhtG)
dhtd = nx.degree_centrality(dhtG)

dhtB = pd.DataFrame.from_dict(dhtb, orient="index", columns=["Betweenness Centrality"])

dhtD = pd.DataFrame.from_dict(dhtd, orient="index", columns=["Degree Centrality"])

dhtB["Degree Centrality"] = dhtD["Degree Centrality"]

nodeClasses = []
for node in dhtB.index:
    nodeClasses += [dhtG.nodes()[node]["nodeClass"]]


dhtB["nodeClass"] = nodeClasses
dhtB["Condition"] = "DHT"


B = pd.concat([ethB, dhtB])


Arbs = B[~((B["nodeClass"] == "upP") | (B["nodeClass"] == "dwP"))]
Arbs = Arbs.sort_values(by=["nodeClass"])

conditions = ["EtOH", "DHT"]
nodeClasses = ["con", "ind", "non"]

boxPairs = [tuple((n,c) for c in conditions) for n in nodeClasses]


len(dhtB)
len(ethB)

fig = plt.figure(figsize=[12,6])

gs = gridspec.GridSpec(ncols=2, nrows=1)

fig.add_subplot(gs[0])

dy="Degree Centrality"; dx="nodeClass"; hue="Condition"
# , "#ACACAC"
colorPalette = ["#5A5A5A", "#F9746D"]

ax=sns.violinplot(x = dx, y = dy, data = Arbs, hue=hue, palette = colorPalette, split=True, inner=None)
ax.legend()
ax1=sns.boxplot( x = dx, y = dy, data = Arbs, hue=hue, color = "black", width = .15, zorder = 10, showcaps = False, boxprops = {'facecolor':'#FFFFFF', "zorder":10}, showfliers=False, whiskerprops = {'linewidth':2, "zorder":10}, saturation = 1, dodge=5)
plt.title("Not Connected \n Degree Centrality")
add_stat_annotation(ax,x=dx, y=dy, data=Arbs, hue=hue,  test="Mann-Whitney", box_pairs=boxPairs, loc='inside' ,line_height=0)


fig.add_subplot(gs[1])

dy="Betweenness Centrality"; dx="nodeClass"; hue="Condition"


ax=sns.violinplot(x = dx, y = dy, data = Arbs, hue=hue, palette = colorPalette, inner=None)
ax.legend()
ax1=sns.boxplot( x = dx, y = dy, data = Arbs, hue=hue, color = "black", width = .15, zorder = 10, showcaps = False, boxprops = {'facecolor':'#FFFFFF', "zorder":10}, showfliers=False, whiskerprops = {'linewidth':2, "zorder":10}, saturation = 1)
plt.title("Not Connected \n Betweenness Centrality")
add_stat_annotation(ax,x=dx, y=dy, data=Arbs, hue=hue,  test="Mann-Whitney", box_pairs=boxPairs, loc='inside', line_height=0)

fig.savefig(f"{figureRoot}/NotConnected.Cicero.Centrality.Violin.pdf")
