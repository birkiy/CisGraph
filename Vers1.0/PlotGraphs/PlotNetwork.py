

from Functions.Packages import *


G = pickle.load(open(f"{dataRoot}/tmpData/GraphsGData.p", "rb" ))


#
# print("Started to draw the SubGraph with gIdx %i" % (gIdx))

nodeClasses = ["upP", "dwP", "con", "ind", "non", "oth"]

pos = nx.spring_layout(G, k=0.2, iterations=20)
#
# for node, ps in pos.items():
#     currentNodeClass = G.nodes()[node]["nodeClass"]
#     nodeClassIdx = nodeClasses.index(currentNodeClass)
#     if nodeClassIdx == 0:  # Down promoter
#         ps[0] += 2
#         ps[1] -= 2
#     elif nodeClassIdx == 1:  # Up Promoter
#         ps[0] -= 2
#         ps[1] -= 2
#     elif nodeClassIdx == 2:  # Con
#         ps[0] += 2
#         ps[1] += 2
#     elif nodeClassIdx == 3:  # Ind
#         ps[0] -= 2
#         ps[1] += 2
#     elif nodeClassIdx == 4:  # Non
#         ps[1] += 3.5
#     elif nodeClassIdx == 5:  # Oth
#         ps[1] -= 3.5
#

plt.close("all")

fig = plt.figure(figsize=[25,25])
######### Draw Graph
nodeClasses = ["upP", "dwP",
               "con", "ind", "non", "oth"]
#
# colors = ["#F95355", "#4C6472",
#           "#FADE89", "#57A4B1", "#B0D894"]
#
# # labels = ["LNCaP DHT Down Promoter",
#           "LNCaP DHT Up Promoter",
#           "Constituvely Active",
#           "Inducible",
#           "Inactive",
#           "Not Clinical"]

for i in nodeClasses:
    plt.plot([], [],
             color=colorPalette[i],
             label=i,
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
    _ = nx.draw_networkx_nodes(G.nodes(),
                           pos=pos,
                           node_size=(G.degree[nodeSingle] * 10) ^ 2,
                           alpha=0.9,
                           nodelist=[nodeSingle],
                           node_color=dictNodeSingle["color"],
                           with_labels=False)


for edgeSingle in G.edges():
    dictEdgeSingle = G.edges()[edgeSingle]
    _ = nx.draw_networkx_edges(G,
                           pos=pos,
                           node_size=10,
                           alpha=0.2,
                           width=math.log(float(dictEdgeSingle["weight"]), 2) * 0.4,
                           arrowsize=4,
                           arrowstyle="->",
                           edgelist=[edgeSingle],
                           edge_color=dictEdgeSingle["color"])

plt.axis("off")

fig.savefig(f"{figureRoot}/PlotNetworkGraph.Final.pdf")
