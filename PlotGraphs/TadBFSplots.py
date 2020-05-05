


from BFS.TadBFS import *


# RootedTAD["tad.X"] is tl
for gIdx in RootedTAD.keys():

    tl = RootedTAD[gIdx]


    gl = [
        val
        for sub in [
            list(T.nodes[t]["subG"].nodes())
            for t in tl
        ]
        for val in sub
    ]

    t = T.subgraph(tl).copy()
    g = G.subgraph(gl).copy()


    for _ in gl:
        if G.nodes[_]["tad"] == "":
            g.remove_node(_)

    posT = nx.drawing.shell_layout(t, nlist=Shell[gIdx])
    posG = nx.circular_layout(g)

    for ps in posG:
        ps2 = posT[G.nodes[ps]["tad"]].copy()
        ps2[0] += random.randint(-65000, 65000) / 1000000
        ps2[1] += random.randint(-65000, 65000) / 1000000

        posG[ps] = ps2




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
                               node_size=(g.degree[nodeSingle] * 30 + 1) ** 1.5,
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

    fig.savefig(home + "/Figures/" + gIdx + "BFS.TADs.pdf")


"""
if len(G.nodes()) < 2:
    continue
pos = nx.drawing.shell_layout(G, nlist=Shell[gIdx])
pos2 = {}
for  key, val in zip(pos.keys(),pos.values()):
    valx = val[0] + random.randint(0,100) / 1000
    valy = val[1] + random.randint(0, 100) / 1000
    pos2[key]= [valx,valy]


fig = plt.figure(figsize=(16, 13))
fig.subplots_adjust(hspace=0.2, wspace=0.2)

ax = fig.add_subplot(1,1,1)

for nodeSingle in G.nodes():
    dictNodeSingle = G.nodes()[nodeSingle]
    nx.draw_networkx_nodes(G.nodes(),
                           pos=pos2,
                           node_size=(G.degree[nodeSingle] * 100) ^ 2,
                           alpha=0.7,
                           nodelist=[nodeSingle],
                           node_color=dictNodeSingle["color"],
                           with_labels=False)

for edgeSingle in G.edges():
    dictEdgeSingle = G.edges()[edgeSingle]

    nx.draw_networkx_edges(G,
                           pos=pos2,
                           node_size=10,
                           alpha=0.7,
                           width=math.log(float(dictEdgeSingle["weight"]), 2) * 0.4,
                           arrowsize=4,
                           arrowstyle="->",
                           edgelist=[edgeSingle],
                           edge_color=dictEdgeSingle["color"])

plt.axis("off")

fig.savefig("GenomicInteractions/Figures/Non/" + gIdx + ".pdf")
"""
