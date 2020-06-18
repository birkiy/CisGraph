

from Functions.GIFunctions import *



dhtG = nx.Graph()
fileG = f"{dataRoot}/GIs/GI.G.DHT.txt"

GnC = fromGI(dhtG, fileG, colorPalette)

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

print("You have an G level graph of %i nodes, %i edges, %i components." % (
    len(dhtG.nodes),
    len(dhtG.edges),
    nx.number_connected_components(dhtG)
))


pickle.dump(dhtG,open(f"{dataRoot}/PickleData/GraphsG.DHT.Data.p", "wb" ))



ethG = nx.Graph()
fileG = f"{dataRoot}/GIs/GI.G.EtOH.txt"

GnC = fromGI(ethG, fileG, colorPalette)

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

print("You have an G level graph of %i nodes, %i edges, %i components." % (
    len(ethG.nodes),
    len(ethG.edges),
    nx.number_connected_components(ethG)
))


pickle.dump(ethG,open(f"{dataRoot}/PickleData/GraphsG.EtOH.Data.p", "wb" ))
