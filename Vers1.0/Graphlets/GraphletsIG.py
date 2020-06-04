

from Functions.Packages import *


ethG = pickle.load(open(f"{dataRoot}/tmpData/GraphsG.EtOH.Data.p", "rb" ))
dhtG = pickle.load(open(f"{dataRoot}/tmpData/GraphsG.DHT.Data.p", "rb" ))



ComponentsE = {
    _:[]
    for _ in range(nx.number_connected_components(ethG))
}
gIdx = 0
tmpOldComp = []
for node in ethG.nodes():
    tmpNewComp = set(nx.node_connected_component(ethG, node))
    if not tmpNewComp in tmpOldComp:
        tmpOldComp.append(tmpNewComp)
        ComponentsE[gIdx] = tmpNewComp
        gIdx += 1

ComponentsD = {
    _:[]
    for _ in range(nx.number_connected_components(dhtG))
}
gIdx = 0
tmpOldComp = []
for node in dhtG.nodes():
    tmpNewComp = set(nx.node_connected_component(dhtG, node))
    if not tmpNewComp in tmpOldComp:
        tmpOldComp.append(tmpNewComp)
        ComponentsD[gIdx] = tmpNewComp
        gIdx += 1

def select_k(spectrum, minimum_energy = 0.9):
    running_total = 0.0
    total = sum(spectrum)
    if total == 0.0:
        return len(spectrum)
    for i in range(len(spectrum)):
        running_total += spectrum[i]
        if running_total / total >= minimum_energy:
            return i + 1
    return len(spectrum)

Similarities = np.empty((len(ComponentsE), len(ComponentsD)), dtype=np.float64)
laplacian1 = []
laplacian2 = []
k1 = []
k2 = []
for keyE in ComponentsE.keys():
    e = ethG.subgraph(ComponentsE[keyE])
    laplacian1.append(nx.spectrum.laplacian_spectrum(e))
    k1.append(select_k(laplacian1[-1]))

for keyD in ComponentsD.keys():
    d = dhtG.subgraph(ComponentsD[keyD])
    laplacian2.append(nx.spectrum.laplacian_spectrum(d))
    k2.append(select_k(laplacian2[-1]))

for keyE in ComponentsE.keys():
    for keyD in ComponentsD.keys():
        k = min(k1[keyE], k2[keyD])
        similarity = sum((laplacian1[keyE][:k] - laplacian2[keyD][:k])**2)
        Similarities[keyE][keyD] = similarity

Similarities = pd.DataFrame(Similarities)


from karateclub import EgoNetSplitter, GraphWave, Role2Vec

eGs = []
for keyE in ComponentsE.keys():
    e = ethG.subgraph(ComponentsE[keyE])
    eGs.append(e)

dGs = []
for keyD in ComponentsD.keys():
    d = dhtG.subgraph(ComponentsD[keyD])
    dGs.append(d)


i = 0
for e in eGs:
    splitter = EgoNetSplitter(1.0)
    mapping = {node: i for node, i in zip(e.nodes(), range(len(e)))}
    e = nx.relabel_nodes(e, mapping)
    splitter.fit(e)
    i += 1
    if splitter.get_memberships()[0] == [] or len(e) == 1:
        continue
    print(splitter.get_memberships(), i-1, len(e))



mapping = {node: i for node, i in zip(e.nodes(), range(len(ethG)))}
eM = nx.relabel_nodes(ethG, mapping)

GW.fit(eM)


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





i = 0
for e in eGs:
    if len(e) == 1:
        continue
    RV = Role2Vec()
    mapping = {node: i for node, i in zip(e.nodes(), range(len(e)))}
    e = nx.relabel_nodes(e, mapping)
    RV.fit(e)
    i += 1
    print(RV.get_embedding(), i-1, len(e))
