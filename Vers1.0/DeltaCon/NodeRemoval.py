

from Function.Packages import *

ethG = pickle.load(open(f"{dataRoot}/tmpData/GraphsG.EtOH.Data.p", "rb" ))
dhtG = pickle.load(open(f"{dataRoot}/tmpData/GraphsG.DHT.Data.p", "rb" ))

dNodes = list(dhtG.nodes())



i = 0
nodeL = []
nodeClassL = []
connectedComponentsL = []
for node in dNodes:
    # df = pd.DataFrame()
    G = dhtG.copy()
    G.remove_node(node)
    nodeL.append(node)
    nodeClassL.append(dhtG.nodes[node]["nodeClass"])
    connectedComponentsL.append(nx.number_connected_components(G))
    # DF = pd.concat([DF,df])

DF = pd.DataFrame()

DF["node"] = nodeL
DF["nodeClass"] = nodeClassL
DF["Connected Components"] = connectedComponentsL

dhtDF = DF



dNodes = list(ethG.nodes())

i = 0
nodeL = []
nodeClassL = []
connectedComponentsL = []
for node in dNodes:
    # df = pd.DataFrame()
    G = ethG.copy()
    G.remove_node(node)
    nodeL.append(node)
    nodeClassL.append(ethG.nodes[node]["nodeClass"])
    connectedComponentsL.append(nx.number_connected_components(G))
    # DF = pd.concat([DF,df])


DF = pd.DataFrame()

DF["node"] = nodeL
DF["nodeClass"] = nodeClassL
DF["Connected Components"] = connectedComponentsL


ethDF = DF



boxPairs = (("con", "ind"),
            ("con", "non"),
            ("con", "dwP"),
            ("con", "upP"),
            ("con", "oth"),
            ("ind", "non"),
            ("ind", "upP"),
            ("ind", "dwP"),
            ("ind", "oth"),
            ("non", "upP"),
            ("non", "oth"),
            ("non", "dwP"),
            ("upP", "dwP"),
            ("upP", "oth"),
            ("oth", "dwP"))


plt.close("all")
fig = plt.figure(figsize=(15,9))
gs = gridspec.GridSpec(ncols=2, nrows=1)

axDht = fig.add_subplot(gs[0])
axDht = sns.violinplot(x="nodeClass", y="Connected Components", data=dhtDF, palette=colorPalette)
# add_stat_annotation(axDht,x="nodeClass", y="Connected Components", data=dhtDF, test="Mann-Whitney", box_pairs=boxPairs, loc='outside', line_height=0)
axDht.set_ylim([2808, 2818])
axEth = fig.add_subplot(gs[1])
axEth = sns.violinplot(x="nodeClass", y="Connected Components", data=ethDF, palette=colorPalette)
# add_stat_annotation(axEth,x="nodeClass", y="Connected Components", data=ethDF, test="Mann-Whitney", box_pairs=boxPairs, loc='outside', line_height=0)
axEth.set_ylim([1930, 1940])
fig.savefig(f"{figureRoot}/nodeRemoval.pdf")
