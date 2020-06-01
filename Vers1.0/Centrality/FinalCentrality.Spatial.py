


from Functions.Packages import *


G = pickle.load(open(f"{dataRoot}/tmpData/GraphsGData.p", "rb" ))



# Not Connected
b = nx.betweenness_centrality(G)
d = nx.degree_centrality(G)

B = pd.DataFrame.from_dict(b, orient="index", columns=["Betweenness Centrality"])

D = pd.DataFrame.from_dict(d, orient="index", columns=["Degree Centrality"])

B["Degree Centrality"] = D["Degree Centrality"]

nodeClasses = []
for node in B.index:
    nodeClasses += [G.nodes()[node]["nodeClass"]]


B["nodeClass"] = nodeClasses

Arbs = B[~((B["nodeClass"] == "upP") | (B["nodeClass"] == "dwP"))]
Arbs = Arbs.sort_values(by=["nodeClass"])

boxPairs = [("con", "ind"),
            ("con", "non"),
            ("ind", "non")]



fig = plt.figure(figsize=[12,6])

gs = gridspec.GridSpec(ncols=2, nrows=1)

fig.add_subplot(gs[0])

dy="Degree Centrality"; dx="nodeClass"; ort="h"; pal = sns.color_palette(n_colors=1)
colorPalette = ["#5A5A5A", "#F9746D", "#ACACAC", "#000000"]

ax=sns.violinplot(x = dx, y = dy, data = Arbs, palette = colorPalette, inner = None)
ax=sns.boxplot( x = dx, y = dy, data = Arbs, color = "black", width = .15, zorder = 10, showcaps = False, boxprops = {'facecolor':'#FFFFFF', "zorder":10}, showfliers=False, whiskerprops = {'linewidth':2, "zorder":10}, saturation = 1)
plt.title("Not Connected \n Degree Centrality")
add_stat_annotation(ax,x=dx, y=dy, data=Arbs,  test="Mann-Whitney", box_pairs=boxPairs, loc='inside' ,line_height=0)



fig.add_subplot(gs[1])

dy="Betweenness Centrality"; dx="nodeClass"; ort="h"; pal = sns.color_palette(n_colors=1)
colorPalette = ["#5A5A5A", "#F9746D", "#ACACAC", "#000000"]

ax=sns.violinplot(x = dx, y = dy, data = Arbs, palette = colorPalette, inner = None)
ax=sns.boxplot( x = dx, y = dy, data = Arbs, color = "black", width = .15, zorder = 10, showcaps = False, boxprops = {'facecolor':'#FFFFFF', "zorder":10}, showfliers=False, whiskerprops = {'linewidth':2, "zorder":10}, saturation = 1)
plt.title("Not Connected \n Betweenness Centrality")
add_stat_annotation(ax,x=dx, y=dy, data=Arbs,  test="Mann-Whitney", box_pairs=boxPairs, loc='inside', line_height=0)

fig.savefig(f"{figureRoot}/NotConnected.Centrality.Violin.pdf")


# Largest Connected
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

ComponentsN = [len(g) for g in ComponentsG.values()]

# Largest => index 4
g = G.subgraph(ComponentsG[4])

b = nx.betweenness_centrality(g)
d = nx.degree_centrality(g)

B = pd.DataFrame.from_dict(b, orient="index", columns=["Betweenness Centrality"])

D = pd.DataFrame.from_dict(d, orient="index", columns=["Degree Centrality"])

B["Degree Centrality"] = D["Degree Centrality"]
nodeClasses = []
for node in B.index:
    nodeClasses += [g.nodes()[node]["nodeClass"]]

B["nodeClass"] = nodeClasses
Arbs = B[~((B["nodeClass"] == "upP") | (B["nodeClass"] == "dwP"))]
Arbs = Arbs.sort_values(by=["nodeClass"])
###########################################


fig = plt.figure(figsize=[12,6])

gs = gridspec.GridSpec(ncols=2, nrows=1)

fig.add_subplot(gs[0])

dy="Degree Centrality"; dx="nodeClass"; ort="h"; pal = sns.color_palette(n_colors=1)
colorPalette = ["#5A5A5A", "#F9746D", "#ACACAC", "#000000"]

ax=sns.violinplot(x = dx, y = dy, data = Arbs, palette = colorPalette, inner = None)
ax=sns.boxplot( x = dx, y = dy, data = Arbs, color = "black", width = .15, zorder = 10, showcaps = False, boxprops = {'facecolor':'#FFFFFF', "zorder":10}, showfliers=False, whiskerprops = {'linewidth':2, "zorder":10}, saturation = 1)
plt.title("Connected \n Degree Centrality")
add_stat_annotation(ax,x=dx, y=dy, data=Arbs,  test="Mann-Whitney", box_pairs=boxPairs, loc='inside' ,line_height=0)



fig.add_subplot(gs[1])

dy="Betweenness Centrality"; dx="nodeClass"; ort="h"; pal = sns.color_palette(n_colors=1)
colorPalette = ["#5A5A5A", "#F9746D", "#ACACAC", "#000000"]

ax=sns.violinplot(x = dx, y = dy, data = Arbs, palette = colorPalette, inner = None)
ax=sns.boxplot( x = dx, y = dy, data = Arbs, color = "black", width = .15, zorder = 10, showcaps = False, boxprops = {'facecolor':'#FFFFFF', "zorder":10}, showfliers=False, whiskerprops = {'linewidth':2, "zorder":10}, saturation = 1)
plt.title("Connected \n Betweenness Centrality")
add_stat_annotation(ax,x=dx, y=dy, data=Arbs,  test="Mann-Whitney", box_pairs=boxPairs, loc='inside', line_height=0)


fig.savefig(f"{figureRoot}/Connected.Centrality.Violin.pdf")
