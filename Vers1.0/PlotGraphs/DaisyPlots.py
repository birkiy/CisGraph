

import MotifSearch.DaisyChains




figH = len(allPaths2) * 3
print(figH)
iInc = 3

fig = plt.figure()
fig = plt.figure(figsize=(20,figH + 2*iInc))
fig.subplots_adjust(hspace=0, wspace=0.2)


# gs = gridspec.GridSpec(1, 4, width_ratios=[1,3,1,1])
# gs = gridspec.GridSpec(1, len(allPaths2))
# print(gs)
# plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0,
#             hspace = iInc, wspace = iInc)
#

i = 0
for path in allPaths2:

    f = int(i / iInc) + 1
    print(f)
    f_ax = fig.add_subplot(len(allPaths2), 1, f)
    #
    # print(path)
    # print("\n")
    # for n in path:
    #
    #     print(list(g.nodes[n]["subP"].nodes()), n)
    gP = g.subgraph(path).copy()

    print("\n")
    #print(gP.nodes(data=True))
    #
    # mapping = {j: path[j] for j in range(len(path))}
    # gP = nx.relabel_nodes(gP, mapping)

    # pl = [
    #     val
    #     for sub in [
    #         list(g.nodes[n]["subP"].nodes())
    #         for n in path
    #     ]
    #     for val in sub
    # ]
    # print(pl)
    # print("\n")
    # p = P.subgraph(pl).copy()
    #

    posG = {path[j] : [-j/3, 0] for j in range(len(path))}
    print(posG)
    print("\n")
    # posP = nx.random_layout(p)
    #
    # for ps in posP:
    #     print(ps, p.nodes[ps]["elm"])
    #     ps2 = posG[p.nodes[ps]["elm"]].copy()
    #     ps2[0] += random.randint(-10000, 10000) / 1000000
    #     ps2[1] += random.randint(-10000, 10000) / 1000000
    #     posP[ps] = ps2



    print("\n")
    """
    pathTypeG = [gP.nodes[_]["nodeClass"] for _ in gP.nodes()]
    pathTypeP = [p.nodes[_]["nodeClass"] for _ in p.nodes()]

    orderG = [colorPalette[x] for x in pathTypeG]
    orderP = [colorPalette[x] for x in pathTypeP]

    caracG = pd.DataFrame({'ID': range(len(pathTypeG)),
                          'myvalue': orderG})
    caracG = carac.set_index("ID")


    caracP = pd.DataFrame({'ID': range(len(pathTypeP)),
                          'myvalue': orderP})
    caracP = carac.set_index("ID")

    """


    for nodeSingle in gP.nodes():
        dictNodeSingle = G.nodes()[nodeSingle]
        nx.draw_networkx_nodes(gP.nodes(),
                               pos=posG,
                               node_size=(len(G.nodes[nodeSingle]["subP"].nodes()) * 40 + G.degree[nodeSingle] * 40) ,
                               alpha=0.8,
                               nodelist=[nodeSingle],
                               node_color=dictNodeSingle["color"],
                               with_labels=False)
    nx.draw_networkx_labels(gP, posG, font_color="#444444")

    posL = {}
    ax = plt.gca()
    for edgeSingle in gP.edges():
        posL[edgeSingle] = []

        edgeLabels = {(u,v): str(round(int(d['distance']) / 1000, 2)) + "kb" for u,v,d in gP.edges(data=True)}

        dictEdgeSingle = G.edges()[edgeSingle]



        pLx = abs(posG[edgeSingle[0]][0] - posG[edgeSingle[1]][0]) / 2
        print(pLx, "x")

        pLy = pLx / 60
        # pLy = (pLx * 1/len(path)**2 * (figH/len(allPaths2)))**2
        print(pLy , "y")

        if posG[edgeSingle[0]][0] - posG[edgeSingle[1]][0] > 0:
            d = posL[edgeSingle] = [posG[edgeSingle[0]][0] - pLx, posG[edgeSingle[0]][1] + pLy + pLx/200]
        else:
            d = posL[edgeSingle] = [posG[edgeSingle[1]][0] - pLx, posG[edgeSingle[1]][1] - pLy - pLx/200]


        l = ax.annotate("",
                xy=posG[edgeSingle[0]], xycoords='data',
                xytext=posG[edgeSingle[1]], textcoords='data',
                arrowprops=dict(arrowstyle="<-", color=dictEdgeSingle["color"],
                                shrinkA=5, shrinkB=5,
                                patchA=None, patchB=None,
                                connectionstyle="arc3,rad=-0.3",
                                ),
                )
        ax.text(posL[edgeSingle][0], posL[edgeSingle][1], edgeLabels[edgeSingle], color="#888888", va="center", ha="center")

    print(posG)
    print(posL)
    #nx.draw_networkx_edge_labels(gP,posL,edge_labels=edge_labels)


    # for nodeSingle in p.nodes():
    #     dictNodeSingle = p.nodes()[nodeSingle]
    #     nx.draw_networkx_nodes(p.nodes(),
    #                             pos=posP,
    #                             node_size=10,
    #                             alpha=0.8,
    #                             nodelist=[nodeSingle],
    #                             node_color=dictNodeSingle["color"],
    #                             with_labels=False)

    """
    nx.draw(gP, posG,
            node_size=300,
            alpha=0.7,
            node_color=caracG['myvalue'],
            edge_alpha=0.5
            )
    nx.draw_networkx_labels(gP, posL, font_size=10)
    nx.draw(p, posP,
            node_size=10,
            )
    """


    i += iInc

    plt.xlim((-0.33 * (len(path)) +0.2, 0.2))

    plt.axis("off")

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
