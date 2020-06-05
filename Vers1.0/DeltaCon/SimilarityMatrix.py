


from Functions.DeltaCon import *






ethG = pickle.load(open(f"{dataRoot}/tmpData/GraphsG.EtOH.Data.p", "rb" ))
dhtG = pickle.load(open(f"{dataRoot}/tmpData/GraphsG.DHT.Data.p", "rb" ))

d, w, E, G1, G2= NX2deltaCon(ethG, dhtG, returnG=True, e=10**-8)


G1 = ethG
G2 = dhtG
nodesU = set(G1.nodes()).union(set(G2.nodes()))

for node in nodesU:
    if not node in G1.nodes():
        att = G2.nodes[node]
        G1.add_node(node)
        for at in att:
            G1.nodes[node][at] = att[at]
    if not node in G2.nodes():
        att = G1.nodes[node]
        G2.add_node(node)
        for at in att:
            G2.nodes[node][at] = att[at]

_G1 = nx.Graph()
_G2 = nx.Graph()
for node in nodesU:
    _G1.add_node(node, nodeClass=G1.nodes[node]["nodeClass"], w= None)
    _G2.add_node(node, nodeClass=G2.nodes[node]["nodeClass"], w= None)

_G1.add_edges_from(list(G1.edges()))
_G2.add_edges_from(list(G2.edges()))

A1 = nx.to_numpy_matrix(_G1)
A2 = nx.to_numpy_matrix(_G2)


d = deltaCon0(A1, A2)
print(f"DeltaCon Similarity score of two graphs is {sim(d)}")

w = deltaConAttrNode(A1, A2)
print(f"Top 20 node impact {sorted(w)[:20]}")
a = sorted([(idx, wi) for idx, wi in zip(range(len(w)), w)], key=lambda kv: kv[1])
b = [_[0] for _ in a]
a = np.array([b]*len(b))
#
A1 = np.array(list(map(lambda x, y: y[x], np.argsort(a), np.array(A1))))
# np.array(list(map(lambda x, y: y[x], np.argsort(b), x)))
A2 = np.array(list(map(lambda x, y: y[x], np.argsort(a), np.array(A2))))




Impact = pd.DataFrame()
Impact["w"] = dict(G1.nodes(data="w")).values()
Impact["nodeClass"] = dict(G1.nodes(data="nodeClass")).values()
Impact["node"] = dict(G1.nodes(data="nodeClass")).keys()



colorPalette = {"con": "#5A5A5A",
                "ind": "#F9746D",
                "non": "#ACACAC",
                "upP": "#000000",
                "dwP": "#000000"}

boxPairs = (("con", "ind"),
            ("con", "non"),
            ("con", "dwP"),
            ("con", "upP"),
            ("ind", "non"),
            ("ind", "upP"),
            ("ind", "dwP"),
            ("non", "upP"),
            ("non", "dwP"),
            ("upP", "dwP"))

fig = plt.figure(figsize=(9,12))
ax = sns.violinplot(x="nodeClass", y="w", data=Impact, palette=colorPalette)
add_stat_annotation(ax,x="nodeClass", y="w", data=Impact,  test="Mann-Whitney", box_pairs=boxPairs, loc='inside', line_height=0)
fig.savefig(f"{figureRoot}/deltaCon.pdf")


plt.close("all")
