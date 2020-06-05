

from karateclub import EgoNetSplitter, GraphWave, Role2Vec
from sklearn.decomposition import PCA


eGs = []
for keyE in ComponentsE.keys():
    e = ethG.subgraph(ComponentsE[keyE])
    eGs.append(e)

dGs = []
for keyD in ComponentsD.keys():
    d = dhtG.subgraph(ComponentsD[keyD])
    dGs.append(d)


EPCA = pd.DataFrame()
DPCA = pd.DataFrame()
for root in both:
    dl = dhtRooted[root]
    d = dhtG.subgraph(dl).copy()
    el = ethRooted[root]
    e = ethG.subgraph(el).copy()
    if len(e) == 1 or len(d) == 1:
        continue
    dGW = GraphWave()
    eGW = GraphWave()
    mapping = {node: i for node, i in zip(d.nodes(), range(len(d)))}
    dM = nx.relabel_nodes(d, mapping)
    mapping = {node: i for node, i in zip(e.nodes(), range(len(e)))}
    eM = nx.relabel_nodes(e, mapping)
    dGW.fit(dM)
    eGW.fit(eM)
    #
    pD = PCA(n_components=2)
    dEmbed = dGW.get_embedding()
    pcaD = pD.fit_transform(dEmbed)
    #
    pE = PCA(n_components=2)
    eEmbed = eGW.get_embedding()
    pcaE = pE.fit_transform(eEmbed)
    #
    # fig = plt.figure(figsize=[7,4])
    # gs = gridspec.GridSpec(ncols=2, nrows=1)
    #
    DpcaDF = pd.DataFrame(data = pcaD, columns = ['principal component 1', 'principal component 2'])
    DpcaDF["nodes"] = list(dict(d.nodes(data="nodeClass")).keys())
    DpcaDF["nodeClass"] = list(dict(d.nodes(data="nodeClass")).values())
    DPCA = pd.concat([DPCA, DpcaDF])
    # ; palette=dColor,
    # dColor = dict(dM.nodes(data="color"))
    # print(dColor)
    # axDht = fig.add_subplot(gs[0])
    # sns.scatterplot(x='principal component 1', y='principal component 2', hue="nodes", palette=dColor, data=DpcaDF, ax=axDht)
    EpcaDF = pd.DataFrame(data = pcaE, columns = ['principal component 1', 'principal component 2'])
    EpcaDF["nodes"] = list(dict(e.nodes(data="nodeClass")).keys())
    EpcaDF["nodeClass"] = list(dict(e.nodes(data="nodeClass")).values())
    EPCA = pd.concat([EPCA, EpcaDF])
    # axEth = fig.add_subplot(gs[1])
    #;  palette=eColor,
    # eColor = dict(eM.nodes(data="color"))
    # print(eColor)
    # sns.scatterplot(x='principal component 1', y='principal component 2', hue="nodes", palette=eColor, data=EpcaDF, ax=axEth)
    #
    # fig.savefig(f"{figureRoot}/CiceroStructurePCA/{root}.pdf")
    # plt.close("all")

EPCA = EPCA.reset_index()
DPCA = DPCA.reset_index()

fig = plt.figure(figsize=[9,6])
gs = gridspec.GridSpec(ncols=2, nrows=1)
axDht = fig.add_subplot(gs[0])
colorPalette = {"ind": "#F9746D", "dwP": "#000000", "upP": "#000000", "non": "#ACACAC","con": "#5A5A5A"}
sns.scatterplot(x='principal component 1', y='principal component 2', hue="nodeClass", palette=colorPalette, data=DPCA, ax=axDht)
axEth = fig.add_subplot(gs[1])
sns.scatterplot(x='principal component 1', y='principal component 2', hue="nodeClass", palette=colorPalette, data=EPCA, ax=axEth)

fig.savefig(f"{figureRoot}/CiceroStructurePCA/total.pdf")


dGW = GraphWave()
eGW = GraphWave()
mapping = {node: i for node, i in zip(dhtG.nodes(), range(len(dhtG)))}
dM = nx.relabel_nodes(dhtG, mapping)
mapping = {node: i for node, i in zip(ethG.nodes(), range(len(ethG)))}
eM = nx.relabel_nodes(ethG, mapping)
dGW.fit(dM)
eGW.fit(eM)
#
pD = PCA(n_components=2)
dEmbed = dGW.get_embedding()
pcaD = pD.fit_transform(dEmbed)
#
pE = PCA(n_components=2)
eEmbed = eGW.get_embedding()
pcaE = pE.fit_transform(eEmbed)






nx.adjacency_matrix(d).todense()





i = 0
for e in eGs:
    if len(e) == 1:
        continue
    GW = GraphWave()
    mapping = {node: i for node, i in zip(e.nodes(), range(len(e)))}
    e = nx.relabel_nodes(e, mapping)
    GW.fit(e)
    i += 1
    # if GW.get_embedding()[0] == [] or len(e) == 1:
    #     continue
    print(GW.get_embedding(), i-1, len(e))
