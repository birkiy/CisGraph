
from Functions.CustomBFS import *




G = dhtG
tTree = {}

Shell = {}
RootedTAD = {}


root="FKBP5"
depthLimit = 4
shell = [None for _ in range(depthLimit+1)]
currentDepth = []

g = bfs_treeCustom(G, root, depth_limit=depthLimit, current_depth=currentDepth, distance=False,shell=shell)


lvlN = depthLimit - min(currentDepth)
if not lvlN in tTree:
    tTree[lvlN] = {}
    if not root in tTree[lvlN]:
        tTree[lvlN][root] = [g]
    else:
        tTree[lvlN][root] += [g]
else:
    if not root in tTree[lvlN]:
        tTree[lvlN][root] = [g]
    else:
        tTree[lvlN][root] += [g]



shell = [i for i in shell if i]


Shell[root] = shell
RootedTAD[root] = list(g.nodes())




gl = RootedTAD[root]

g = G.subgraph(gl).copy()

posG = nx.drawing.shell_layout(g, nlist=Shell[root])

[5 for _ in g.nodes(data="nodeClass") if _[1] == "upP"]





plt.close("all")

fig = plt.figure(figsize=(15,15))
fig.subplots_adjust(hspace=0.2, wspace=0.2)
ax = fig.add_subplot(1, 1, 1)



for nodeSingle in g.nodes():
    dictNodeSingle = g.nodes()[nodeSingle]
    nx.draw_networkx_nodes(g.nodes(),
                           pos=posG,
                           node_size=(G.degree[nodeSingle] * 100) ,
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


fig.savefig(f"{figureRoot}/{root}.5lvl.DHT.BFS.pdf")
