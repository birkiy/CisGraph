

print("Started to draw the SubGraph with gIdx %i" % (gIdx))

nodeClasses = ["upP", "dwP",
               "con", "ind", "non"]

pos = nx.spring_layout(G, k=0.2, iterations=20)

for node, ps in pos.items():
    currentNodeClass = G.nodes()[node]["nodeClass"]
    if currentNodeClass in nodeClasses:
        nodeClassIdx = nodeClasses.index(currentNodeClass)
    if nodeClassIdx == 0:  # Down promoter
        ps[0] += 2.5
        ps[1] -= 3
    elif nodeClassIdx == 1:  # Up Promoter
        ps[0] -= 2.5
        ps[1] -= 3
    elif nodeClassIdx == 2:  # Con
        ps[0] += 4.5
        ps[1] += 1
    elif nodeClassIdx == 3:  # Ind
        ps[0] -= 4.5
        ps[1] += 1
    elif nodeClassIdx == 4:  # Non
        ps[1] += 4.5




fig = plt.figure(figsize=[25,20])
######### Draw Graph
locs = ["lower left", "lower right", "center right", "center left", "upper center"]

nodeClasses = ["upP", "dwP",
               "con", "ind", "non"]

colors = ["#F95355", "#4C6472",
          "#FADE89", "#57A4B1", "#B0D894"]

labels = ["LNCaP DHT Down Promoter",
          "LNCaP DHT Up Promoter",
          "Constituvely Active",
          "Inducible",
          "Inactive"]

for i in range(0, len(colors)):
    plt.plot([], [],
             color=colors[i],
             label=labels[i],
             alpha=0.5,
             marker=".",
             markersize=25
             )


plt.legend(loc="lower center",
           frameon=False,
           scatterpoints=1,
           mode="expand",
           bbox_to_anchor=(-0.1, 0.8),
           fontsize=30
           )

for nodeSingle in G.nodes():
    dictNodeSingle = G.nodes()[nodeSingle]
    nx.draw_networkx_nodes(G.nodes(),
                           pos=pos,
                           node_size=(G.degree[nodeSingle] * 10) ^ 2,
                           alpha=0.7,
                           nodelist=[nodeSingle],
                           node_color=dictNodeSingle["color"],
                           with_labels=False)


for edgeSingle in G.edges():
    dictEdgeSingle = G.edges()[edgeSingle]
    nx.draw_networkx_edges(G,
                           pos=pos,
                           node_size=10,
                           alpha=0.2,
                           width=math.log(float(dictEdgeSingle["weight"]), 2) * 0.4,
                           arrowsize=4,
                           arrowstyle="->",
                           edgelist=[edgeSingle],
                           edge_color=dictEdgeSingle["color"])

plt.axis("off")

fig.savefig("PlotNetworkGraph.pdf")
