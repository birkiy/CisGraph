

from Functions.GIFunctions import *



G = nx.Graph()
fileG = f"{dataRoot}/GIs/GI.G.txt"

GnC = fromGI(G, fileG, colorPalette)

ComponentsG = {
    _:[]
    for _ in range(nx.number_connected_components(G))
}
gIdx = 0
tmpOldComp = []
for node in G.nodes():
    tmpNewComp = set(nx.node_connected_component(G, node))
    if not tmpNewComp in tmpOldComp:
        tmpOldComp.append(tmpNewComp)
        ComponentsG[gIdx] = tmpNewComp
        gIdx += 1

print("You have an G level graph of %i nodes, %i edges, %i components." % (
    len(G.nodes),
    len(G.edges),
    nx.number_connected_components(G)
))




T = nx.Graph()
fileT = f"{dataRoot}/GIs/GI.T.txt"

TnC = fromGI(T, fileT, colorPalette)

ComponentsT = {
    _:[]
    for _ in range(nx.number_connected_components(T))
}
gIdx = 0
tmpOldComp = []
for node in T.nodes():
    tmpNewComp = set(nx.node_connected_component(T, node))
    if not tmpNewComp in tmpOldComp:
        tmpOldComp.append(tmpNewComp)
        ComponentsT[gIdx] = tmpNewComp
        gIdx += 1

print("You have an T level graph of %i nodes, %i edges, %i components." % (
    len(T.nodes),
    len(T.edges),
    nx.number_connected_components(T)
))





C = nx.Graph()
fileC = f"{dataRoot}/GIs/GI.C.txt"

CnC = fromGI(C, fileC, colorPalette)

ComponentsC = {
    _:[]
    for _ in range(nx.number_connected_components(C))
}
gIdx = 0
tmpOldComp = []
for node in C.nodes():
    tmpNewComp = set(nx.node_connected_component(C, node))
    if not tmpNewComp in tmpOldComp:
        tmpOldComp.append(tmpNewComp)
        ComponentsC[gIdx] = tmpNewComp
        gIdx += 1

print("You have an C level graph of %i nodes, %i edges, %i components." % (
    len(C.nodes),
    len(C.edges),
    nx.number_connected_components(C)
))



pickle.dump(ComponentsC,open(f"{dataRoot}/PickleData/ComponentsCData.p", "wb" ))
pickle.dump(ComponentsT,open(f"{dataRoot}/PickleData/ComponentsTData.p", "wb" ))
pickle.dump(ComponentsG,open(f"{dataRoot}/PickleData/ComponentsGData.p", "wb" ))

pickle.dump(C,open(f"{dataRoot}/PickleData/GraphsCData.p", "wb" ))
pickle.dump(T,open(f"{dataRoot}/PickleData/GraphsTData.p", "wb" ))
pickle.dump(G,open(f"{dataRoot}/PickleData/GraphsGData.p", "wb" ))
