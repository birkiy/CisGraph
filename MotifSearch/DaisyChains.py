
from Functions.Helpers import *

from Functions.Packages import *

home = "/home/birkiy/github/CisGraph"

C = pickle.load(open(home + "/Data/tmpData/GraphsCData.p", "rb" ))
T = pickle.load(open(home + "/Data/tmpData/GraphsTData.p", "rb" ))
G = pickle.load(open(home + "/Data/tmpData/GraphsGData.p", "rb" ))



colorPalette ={"hom": "#888888",
               "con": "#FADE89",
               "ind": "#57A4B1",
               "non": "#B0D894",
               "tad": "#4c5172",
               "chr": "#5a4c72",
               "upP": "#F95355",
               "dwP": "#4C6472"
}




upP = [_[0] for _ in G.nodes(data="nodeClass") if _[1] == "upP"]






source = "FKBP5"

rangeX = G.nodes[source]["nodeRange"]
print(range)
Lbed = dict(G.nodes(data="nodeRange"))

Lbed = {k: v for k,v in sorted(Lbed.items(), key=lambda kv: kv[1][1])}

n = rangesFromUpperRange(rangeX[0], rangeX[1] - 200*1000, rangeX[2] + 200*1000,[Lbed],["L"] )

n = [_[0] for _ in n if not G.nodes[_[0]]["nodeClass"] in ("hom", "upP", "dwP")]

n.append(source)

g = G.subgraph(n).copy()


allPaths = []
for node in g.nodes():

    paths = list(nx.all_simple_paths(g,source, node))
    allPaths += paths



fig = plt.figure(figsize=(30,20))
# gs = gridspec.GridSpec(1, 4, width_ratios=[1,3,1,1])
gs = gridspec.GridSpec(1, 1, width_ratios=[1])
ax0 = plt.subplot(gs[0])


i = 0
for path in allPaths:


    p = nx.path_graph(len(path))
    mapping = {j: path[j] for j in range(len(path))}
    print(path)
    print(mapping)
    r = nx.relabel_nodes(p, mapping)

    pos = {path[j] : [-j/3,-i] for j in range(len(path))}
    pos2 = {path[j] : [-j/3,-i+.1] for j in range(len(path))}
    print(pos)

    pathType = [G.nodes[_]["nodeClass"] for _ in path]

    order = [colorPalette[x] for x in pathType]

    carac = pd.DataFrame({'ID': range(len(pathType)),
                          'myvalue': order})
    carac = carac.set_index("ID")

    nx.draw(r, pos,
            node_size=300,
            alpha=0.7,
            node_color=carac['myvalue'],
            edge_alpha=0.5
            )
    nx.draw_networkx_labels(r, pos2,font_size=10)
    i += 1

    plt.xlim((-len(path), 1))

fig.savefig(home + "/Figures/" + source + ".pdf" )

plt.close("all")




"""
print(G.edges[('FKBP5', 'overlapped_read_328-ARBS')]["distance"], G.nodes['overlapped_read_328-ARBS']["nodeClass"])


print(G.edges[('FKBP5', 'peakno_667-ARBS')]["distance"], G.nodes['peakno_667-ARBS']["nodeClass"])
print(G.edges[('FKBP5', 'peakno_672-ARBS')]["distance"], G.nodes['peakno_672-ARBS']["nodeClass"])
print(G.edges[('FKBP5', 'overlapped_read_327-ARBS')]["distance"], G.nodes['overlapped_read_327-ARBS']["nodeClass"])

allPaths = nx.single_source_shortest_path(G,"FKBP5")

for target in allPaths.keys():
    path = allPaths[target]

    condition = [
        1 if
          G.edges[(path[idx-1], path[idx])]["distance"] == "INF"
          or abs(float(G.edges[(path[idx-1], path[idx])]["distance"])) > 100 * 1000
          or G.nodes[path[idx]]["nodeClass"] in ("hom", "upP", "dwP")
        else 0
        for idx in range(1, len(path))
    ]
    try:
        aa = condition.index(1)
        path = path[0:aa+1]
    except:
        pass
    allPaths[target] = path



allPathsList = list(allPaths.values())
print(allPaths)
allPaths = []
for path in allPathsList:
    if not path in allPaths and len(path) > 1:
        allPaths += [path]

print(allPaths)

fig = plt.figure(figsize=(15,20))
gs = gridspec.GridSpec(1, 4, width_ratios=[1,3,1,1])
ax0 = plt.subplot(gs[1])


i = 0
for path in allPaths:

    p = nx.path_graph(len(path))
    pos = {j : [-j/3,-i] for j in range(len(path))}

    pathType = [G.nodes[_]["nodeClass"] for _ in path]

    order = [colorPalette[x] for x in pathType]

    carac = pd.DataFrame({'ID': range(len(pathType)),
                          'myvalue': order})
    carac = carac.set_index("ID")

    nx.draw(p, pos,
            node_size=300,
            alpha=0.7,
            node_color=carac['myvalue'],
            edge_alpha=0.5)
    i += 1

    if i == 50:
        break

plt.show()



ne = G.neighbors("FKBP5")
for n in ne:
    print(G.edges[("FKBP5", n)]["distance"],G.nodes[n]["nodeClass"] , n)

print("")
ne = G.neighbors("overlapped_read_328-ARBS")
for n in ne:
    print(G.edges[("overlapped_read_328-ARBS", n)]["distance"],G.nodes[n]["nodeClass"] , n)


print("")
ne = G.neighbors("peakno_674-ARBS")
for n in ne:
    print(G.edges[("peakno_674-ARBS", n)]["distance"],G.nodes[n]["nodeClass"] , n)



# print(g.nodes(data=True))
print(list(nx.all_simple_paths(g, "FKBP5", "overlapped_read_328-ARBS")))

print("1")
print(list(nx.all_simple_paths(g, "FKBP5", "peakno_674-ARBS")))

print("2")
print(list(nx.shortest_path(G, "FKBP5", "overlapped_read_329-ARBS")))

print("3")
print(g.edges())
"""
