
from Functions.Packages import *


C = pickle.load(open(home + "/Data/tmpData/GraphsCData.p", "rb" ))
T = pickle.load(open(home + "/Data/tmpData/GraphsTData.p", "rb" ))
G = pickle.load(open(home + "/Data/tmpData/GraphsGData.p", "rb" ))




tl = list(T.neighbors(G.nodes["STEAP4"]["tad"]))
print(tl)

gl = [
    val
    for sub in [
        list(T.nodes[t]["subG"].nodes())
        for t in tl
    ]
    for val in sub
]


pl = [
    val
    for sub in [
        list(G.nodes[t]["subP"].nodes())
        for t in gl
    ]
    for val in sub
]


t = T.subgraph(tl).copy()
g = G.subgraph(gl).copy()
p = P.subgraph(pl).copy()


for _ in gl:
    if G.nodes[_]["tad"] == "":
        g.remove_node(_)

posT = nx.fruchterman_reingold_layout(t)
posG = nx.random_layout(g)
posP = nx.random_layout(p)


for ps in posG:
    ps2 = posT[G.nodes[ps]["tad"]].copy()
    ps2[0] += random.randint(-65000, 65000) / 1000000
    ps2[1] += random.randint(-65000, 65000) / 1000000

    posG[ps] = ps2

for ps in posP:
    ps2 = posG[P.nodes[ps]["elm"]].copy()
    ps2[0] += random.randint(-65000, 65000) / 1000000
    ps2[1] += random.randint(-65000, 65000) / 1000000

    posP[ps] = ps2





plt.close("all")

fig = plt.figure(figsize=(50, 50))
fig.subplots_adjust(hspace=0.2, wspace=0.2)
ax = fig.add_subplot(1, 1, 1)




for nodeSingle in t.nodes():
    dictNodeSingle = t.nodes()[nodeSingle]
    nx.draw_networkx_nodes(t.nodes(),
                           pos=posT,
                           node_size=(t.degree[nodeSingle] * 600 + 1) ** 1.5,
                           alpha=0.2,
                           nodelist=[nodeSingle],
                           node_color=dictNodeSingle["color"],
                           with_labels=False)
"""
for edgeSingle in t.edges():
    dictEdgeSingle = t.edges()[edgeSingle]

    nx.draw_networkx_edges(t,
                           pos=posT,
                           node_size=15,
                           alpha=0.5,
                           width=math.log(float(dictEdgeSingle["weight"]), 2) ,
                           arrowsize=4,
                           arrowstyle="->",
                           edgelist=[edgeSingle],
                           edge_color=dictEdgeSingle["color"])
    # dictEdgeSingle["color"]
"""



for nodeSingle in g.nodes():
    dictNodeSingle = g.nodes()[nodeSingle]
    nx.draw_networkx_nodes(g.nodes(),
                           pos=posG,
                           node_size=len(g.nodes[nodeSingle]["subP"].nodes()) * 30 + 1) ** 1.5,
                           alpha=0.8,
                           nodelist=[nodeSingle],
                           node_color=dictNodeSingle["color"],
                           with_labels=False)

for edgeSingle in g.edges():
    dictEdgeSingle = g.edges()[edgeSingle]

    nx.draw_networkx_edges(g,
                           pos=posG,
                           node_size=15,
                           alpha=0.5,
                           width=math.log(float(dictEdgeSingle["weight"]), 2) * 0.4,
                           arrowsize=4,
                           arrowstyle="->",
                           edgelist=[edgeSingle],
                           edge_color=dictEdgeSingle["color"])
    # dictEdgeSingle["color"]



for nodeSingle in p.nodes():
    dictNodeSingle = p.nodes()[nodeSingle]
    nx.draw_networkx_nodes(p.nodes(),
                           pos=posP,
                           node_size=10,
                           alpha=0.8,
                           nodelist=[nodeSingle],
                           node_color=dictNodeSingle["color"],
                           with_labels=False)




plt.axis("off")

fig.savefig(home + "/Figures/" + "KLK3" + ".TAD.pdf")
