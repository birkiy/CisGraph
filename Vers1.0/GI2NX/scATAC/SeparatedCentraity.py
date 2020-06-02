




ethG = pickle.load(open(f"{dataRoot}/tmpData/GraphsG.EtOH.Data.p", "rb" ))
dhtG = pickle.load(open(f"{dataRoot}/tmpData/GraphsG.DHT.Data.p", "rb" ))

Gs = [dhtG, ethG]


G = Gs[0]
###################################3
ComponentsG = {
    _:[]
    for _ in range(nx.number_connected_components(G))
}
gIdx = 0
tmpOldComp = []
for node in G.nodes():
    tmpNewComp = set(nx.node_connected_component(G, node))
    if not tmpNewComp in tmpOldComp:
        tmpOldComp.append(tmpNewComp)
        ComponentsG[gIdx] = tmpNewComp
        gIdx += 1


dhtDF = pd.DataFrame()
for gIdx in ComponentsG.keys():
    g = G.subgraph(ComponentsG[gIdx])
    if len(g) == 1:
        continue
    b = nx.betweenness_centrality(g)
    d = nx.degree_centrality(g)
    B = pd.DataFrame.from_dict(b, orient="index", columns=["Betweenness Centrality"])
    D = pd.DataFrame.from_dict(d, orient="index", columns=["Degree Centrality"])
    B["Degree Centrality"] = D["Degree Centrality"]
    nodeClasses = []
    for node in B.index:
        nodeClasses += [g.nodes()[node]["nodeClass"]]
    B["nodeClass"] = nodeClasses
    B["Condition"] = "DHT"
    dhtDF = pd.concat([dhtDF, B])


#########################################3


G = Gs[1]
###################################3
ComponentsG = {
    _:[]
    for _ in range(nx.number_connected_components(G))
}
gIdx = 0
tmpOldComp = []
for node in G.nodes():
    tmpNewComp = set(nx.node_connected_component(G, node))
    if not tmpNewComp in tmpOldComp:
        tmpOldComp.append(tmpNewComp)
        ComponentsG[gIdx] = tmpNewComp
        gIdx += 1

ethDF = pd.DataFrame()
for gIdx in ComponentsG.keys():
    g = G.subgraph(ComponentsG[gIdx])
    if len(g) == 1:
        continue
    b = nx.betweenness_centrality(g)
    d = nx.degree_centrality(g)
    B = pd.DataFrame.from_dict(b, orient="index", columns=["Betweenness Centrality"])
    D = pd.DataFrame.from_dict(d, orient="index", columns=["Degree Centrality"])
    B["Degree Centrality"] = D["Degree Centrality"]
    nodeClasses = []
    for node in B.index:
        nodeClasses += [g.nodes()[node]["nodeClass"]]
    B["nodeClass"] = nodeClasses
    B["Condition"] = "EtOH"
    ethDF = pd.concat([ethDF, B])


#########################################3


DF = pd.concat([ethDF, dhtDF])


Arbs = DF[~((DF["nodeClass"] == "upP") | (DF["nodeClass"] == "dwP"))]
Arbs = Arbs.sort_values(by=["nodeClass"])








fig = plt.figure(figsize=[12,6])

gs = gridspec.GridSpec(ncols=2, nrows=1)

fig.add_subplot(gs[0])

dy="Degree Centrality"; dx="nodeClass"; hue="Condition"
# , "#ACACAC"
colorPalette = ["#5A5A5A", "#F9746D"]

ax=sns.violinplot(x = dx, y = dy, data = Arbs, hue=hue, palette = colorPalette)
ax.legend()
# ax1=sns.boxplot( x = dx, y = dy, data = Arbs, hue=hue, color = "black", width = .15, zorder = 10, showcaps = False, boxprops = {'facecolor':'#FFFFFF', "zorder":10}, showfliers=False, whiskerprops = {'linewidth':2, "zorder":10}, saturation = 1, dodge=5)
plt.title("Each Connected \n Degree Centrality")
# add_stat_annotation(ax,x=dx, y=dy, data=Arbs, hue=hue,  test="Mann-Whitney", box_pairs=boxPairs, loc='inside' ,line_height=0)


fig.add_subplot(gs[1])

dy="Betweenness Centrality"; dx="nodeClass"; hue="Condition"


ax=sns.violinplot(x = dx, y = dy, data = Arbs, hue=hue, palette = colorPalette)
ax.legend()
# ax1=sns.boxplot( x = dx, y = dy, data = Arbs, hue=hue, color = "black", width = .15, zorder = 10, showcaps = False, boxprops = {'facecolor':'#FFFFFF', "zorder":10}, showfliers=False, whiskerprops = {'linewidth':2, "zorder":10}, saturation = 1)
plt.title("Each Connected \n Betweenness Centrality")
# add_stat_annotation(ax,x=dx, y=dy, data=Arbs, hue=hue,  test="Mann-Whitney", box_pairs=boxPairs, loc='inside', line_height=0)

fig.savefig(f"{figureRoot}/EachConnected.Cicero.Centrality.Violin.pdf")
