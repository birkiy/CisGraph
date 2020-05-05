


from Functions.Packages import *

home = "/home/birkiy/github/CisGraph"

C = pickle.load(open(home + "/Data/tmpData/GraphsCData.p", "rb" ))
T = pickle.load(open(home + "/Data/tmpData/GraphsTData.p", "rb" ))
G = pickle.load(open(home + "/Data/tmpData/GraphsGData.p", "rb" ))

"""
print(G.nodes["STEAP4"]["chr"])

G.nodes["STEAP4"]["chr"]
G.nodes["STEAP4"]["chr"]
"""



t = C
g = T

gl = list(g.nodes()).copy()

for _ in gl:
    if T.nodes[_]["chr"] == "":
        g.remove_node(_)

posT = nx.spring_layout(t, k=0.9, iterations=20)
posG = nx.circular_layout(g)

for ps in posG:
    ps2 = posT[T.nodes[ps]["chr"]].copy()
    ps2[0] += random.randint(-70000, 70000) / 1000000
    ps2[1] += random.randint(-70000, 70000) / 1000000

    posG[ps] = ps2




plt.close("all")

fig = plt.figure(figsize=(50, 50))
fig.subplots_adjust(hspace=0.2, wspace=0.2)
ax = fig.add_subplot(1, 1, 1)




for nodeSingle in t.nodes():
    dictNodeSingle = t.nodes()[nodeSingle]
    nx.draw_networkx_nodes(t.nodes(),
                           pos=posT,
                           node_size=(t.degree[nodeSingle] * 300 + 1) ** 1.5,
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
                           node_size=(g.degree[nodeSingle] * 50),
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

plt.axis("off")

fig.savefig(home + "/Figures/" + "All" + "Chrom.pdf")
