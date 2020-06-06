


from Functions.DeltaCon import *



ethG = pickle.load(open(f"{dataRoot}/tmpData/GraphsG.EtOH.Data.p", "rb" ))
dhtG = pickle.load(open(f"{dataRoot}/tmpData/GraphsG.DHT.Data.p", "rb" ))




ComponentsG = {
    _:[]
    for _ in range(nx.number_connected_components(dhtG))
}
gIdx = 0
tmpOldComp = []
for node in dhtG.nodes():
    tmpNewComp = set(nx.node_connected_component(dhtG, node))
    if not tmpNewComp in tmpOldComp:
        tmpOldComp.append(tmpNewComp)
        ComponentsG[gIdx] = tmpNewComp
        gIdx += 1

ComponentsN = [len(g) for g in ComponentsG.values()]



dhtSing = [ComponentsG[_] for _ in ComponentsG.keys() if len(ComponentsG[_]) == 1 ]


ComponentsG = {
    _:[]
    for _ in range(nx.number_connected_components(ethG))
}
gIdx = 0
tmpOldComp = []
for node in ethG.nodes():
    tmpNewComp = set(nx.node_connected_component(ethG, node))
    if not tmpNewComp in tmpOldComp:
        tmpOldComp.append(tmpNewComp)
        ComponentsG[gIdx] = tmpNewComp
        gIdx += 1

ComponentsN = [len(g) for g in ComponentsG.values()]

ethSing = [ComponentsG[_] for _ in ComponentsG.keys() if len(ComponentsG[_]) == 1 ]





ethSingle = set()
for node in ethSing:
    ethSingle = ethSingle.union(node)



dhtSingle = set()
for node in dhtSing:
    dhtSingle = dhtSingle.union(node)



BothSingle = dhtSingle.intersection(ethSingle)


for node in BothSingle:
    ethG.remove_node(node)
    dhtG.remove_node(node)




d, w, E, G1, G2= NX2deltaCon(ethG, dhtG, returnG=True)



Impact = pd.DataFrame()
Impact["w"] = dict(G1.nodes(data="w")).values()
Impact["nodeClass"] = dict(G1.nodes(data="nodeClass")).values()
Impact["node"] = dict(G1.nodes(data="nodeClass")).keys()


Arbs = Impact[~((Impact["nodeClass"] == "upP") | (Impact["nodeClass"] == "dwP"))]

boxPairs = (("con", "ind"),
            ("con", "non"),
            ("con", "oth"),
            ("ind", "non"),
            ("ind", "oth"),
            ("non", "oth"))


fig = plt.figure(figsize=(5,8))
ax = sns.boxplot(x="nodeClass", y="w", data=Arbs, palette=colorPalette, width = 0.5)
add_stat_annotation(ax,x="nodeClass", y="w", data=Arbs,  test="Mann-Whitney", box_pairs=boxPairs, loc='inside', line_height=0)
plt.title("DHT vs. EtOH Co-accesibility\n (DeltaCon Scores)")
fig.savefig(f"{figureRoot}/deltaConE.R.DHTvsEtOH.white.pdf")


plt.close("all")






ethG

dhtG
