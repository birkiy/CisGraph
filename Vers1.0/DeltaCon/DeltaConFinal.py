


from Functions.DeltaCon import *



ethG = pickle.load(open(f"{dataRoot}/tmpData/GraphsG.EtOH.Data.p", "rb" ))
dhtG = pickle.load(open(f"{dataRoot}/tmpData/GraphsG.DHT.Data.p", "rb" ))



#
# ComponentsG = {
#     _:[]
#     for _ in range(nx.number_connected_components(dhtG))
# }
# gIdx = 0
# tmpOldComp = []
# for node in dhtG.nodes():
#     tmpNewComp = set(nx.node_connected_component(dhtG, node))
#     if not tmpNewComp in tmpOldComp:
#         tmpOldComp.append(tmpNewComp)
#         ComponentsG[gIdx] = tmpNewComp
#         gIdx += 1
#
# ComponentsN = [len(g) for g in ComponentsG.values()]
#
#
#
# dhtSing = [ComponentsG[_] for _ in ComponentsG.keys() if len(ComponentsG[_]) == 1 ]
#
#
# ComponentsG = {
#     _:[]
#     for _ in range(nx.number_connected_components(ethG))
# }
# gIdx = 0
# tmpOldComp = []
# for node in ethG.nodes():
#     tmpNewComp = set(nx.node_connected_component(ethG, node))
#     if not tmpNewComp in tmpOldComp:
#         tmpOldComp.append(tmpNewComp)
#         ComponentsG[gIdx] = tmpNewComp
#         gIdx += 1
#
# ComponentsN = [len(g) for g in ComponentsG.values()]
#
# ethSing = [ComponentsG[_] for _ in ComponentsG.keys() if len(ComponentsG[_]) == 1 ]
#
#
#
#
#
# ethSingle = set()
# for node in ethSing:
#     ethSingle = ethSingle.union(node)
#
#
#
# dhtSingle = set()
# for node in dhtSing:
#     dhtSingle = dhtSingle.union(node)
#
#
#
# BothSingle = dhtSingle.intersection(ethSingle)
#
#
# for node in BothSingle:
#     ethG.remove_node(node)
#     dhtG.remove_node(node)
#
#
#

d, w, E, G1, G2= NX2deltaCon(ethG, dhtG, returnG=True)



Impact = pd.DataFrame()
Impact["w"] = dict(G1.nodes(data="w")).values()
Impact["nodeClass"] = dict(G1.nodes(data="nodeClass")).values()
Impact["node"] = dict(G1.nodes(data="nodeClass")).keys()


Arbs = Impact[~((Impact["nodeClass"] == "upP") | (Impact["nodeClass"] == "dwP") | (Impact["nodeClass"] == "oth"))]

boxPairs = (("con", "ind"),
            ("con", "non"),
            ("ind", "non"))


fig = plt.figure(figsize=(6,8))
ax = sns.violinplot(x="nodeClass", y="w", data=Arbs, palette=colorPalette, width = 0.5, order = ["con", "ind", "non"], flierprops={"marker": "."})
# ax = sns.boxplot( x = "nodeClass", y = "w", data = Arbs, color = "black", width = .05, zorder = 10, showcaps = False, boxprops = {'facecolor':'#FFFFFF', "zorder":10}, showfliers=False, whiskerprops = {'linewidth':2, "zorder":10}, saturation = 1)
add_stat_annotation(ax,x="nodeClass", y="w", data=Arbs,  test="Mann-Whitney", box_pairs=boxPairs, loc='inside', line_height=0)
ax.set_xlabel('Node Classes', fontsize=16)
ax.set_ylabel('w', fontsize=16)
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(14)

plt.title("DHT vs. EtOH Co-accesibility\n (DeltaCon Scores)")
fig.savefig(f"{figureRoot}/deltaConE.DHTvsEtOH.Final.violin.pdf")


plt.close("all")



ImpactE = pd.DataFrame()
nodes = list(G1.nodes())

for eI in E:
    src = eI[0]
    tar = eI[1]
    e1 = G1.nodes[nodes[src]]["nodeClass"]
    e2 = G1.nodes[nodes[tar]]["nodeClass"]
    dir = eI[2]
    score = eI[3]
    if 0:
        di = {"edge": [(src, tar), (tar, src)], "edgeType": [(e1, e2), (e2, e1)], "w": [score, score], "Direction": [dir, dir] }
    else:
        di = {"edge": [(src, tar)], "edgeType": [(e1, e2)], "w": [score], "Direction": [dir] }
    tmp = pd.DataFrame(di)
    ImpactE = pd.concat([ImpactE, tmp])




x="edgeType",

boxPairs = [("+", "-")]
fig = plt.figure(figsize=(30,8))
ax = sns.violinplot(y="w", x="Direction", data=ImpactF, width = 0.5)
add_stat_annotation(ax, y="w", x="Direction", data=ImpactF,  test="Mann-Whitney", box_pairs=boxPairs, loc='inside', line_height=0)
# plt.title("DHT vs. EtOH Co-accesibility\n (DeltaCon Scores)")
plt.xticks(rotation=45)
fig.savefig(f"{figureRoot}/deltaConE.DHTvsEtOH.Edge.pdf")


plt.close("all")









edTy = [('non', 'oth'),
       ('non', 'non'),
       ('non', 'upP'),
       ('non', 'ind'),
       ('non', 'con'),
       ('non', 'dwP'),
       ('oth', 'oth'),
       ('oth', 'con'),
       ('oth', 'upP'),
       ('oth', 'dwP'),
       ('oth', 'ind'),
       ('upP', 'upP'),
       ('upP', 'dwP'),
       ('upP', 'con'),
       ('upP', 'ind'),
       ('con', 'dwP'),
       ('con', 'con'),
       ('con', 'ind'),
       ('dwP', 'dwP'),
       ('dwP', 'ind'),
       ('ind', 'ind')]

ImpactF = ImpactE[ImpactE["edgeType"].isin(edTy)]



boxPairs = [((_, "+"), (_, "-")) for _ in ImpactF["edgeType"].unique()]


fig = plt.figure(figsize=(30,8))
ax = sns.violinplot(x="edgeType", y="w", hue="Direction", data=ImpactF, width = 0.5)
add_stat_annotation(ax,x="edgeType", y="w", hue="Direction", data=ImpactF,  test="Mann-Whitney", box_pairs=boxPairs, loc='inside', line_height=0)
# plt.title("DHT vs. EtOH Co-accesibility\n (DeltaCon Scores)")
plt.xticks(rotation=45)
fig.savefig(f"{figureRoot}/deltaConE.DHTvsEtOH.Edge.pdf")


plt.close("all")


ethG

dhtG
