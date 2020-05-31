

from Functions.GIFunctions import *





ethG = nx.Graph()
fileG = home + "/Data/GIs/GI.G.scATAC.ETOH.txt"

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

print("You have an G EtOH level graph of %i nodes, %i edges, %i components." % (
    len(ethG.nodes),
    len(ethG.edges),
    nx.number_connected_components(ethG)
))

pickle.dump(ethG,open(home + "/Data/tmpData/Graphs.EtOH.GData.p", "wb" ))


dhtG = nx.Graph()
fileG = home + "/Data/GIs/GI.G.scATAC.DHT.txt"

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

print("You have an G DHT level graph of %i nodes, %i edges, %i components." % (
    len(dhtG.nodes),
    len(dhtG.edges),
    nx.number_connected_components(dhtG)
))

pickle.dump(dhtG,open(home + "/Data/tmpData/GraphsG.DHT.Data.p", "wb" ))


# T = nx.Graph()
# fileT = home + "/Data/GIs/GI.T.txt"
#
# TnC = fromGI(T, fileT, colorPalette)
#
# ComponentsT = {
#     _:[]
#     for _ in range(nx.number_connected_components(T))
# }
# gIdx = 0
# tmpOldComp = []
# for node in T.nodes():
#     tmpNewComp = set(nx.node_connected_component(T, node))
#     if not tmpNewComp in tmpOldComp:
#         tmpOldComp.append(tmpNewComp)
#         ComponentsT[gIdx] = tmpNewComp
#         gIdx += 1
#
# print("You have an T level graph of %i nodes, %i edges, %i components." % (
#     len(T.nodes),
#     len(T.edges),
#     nx.number_connected_components(T)
# ))
#
#
#
#
#
# C = nx.Graph()
# fileC = home + "/Data/GIs/GI.C.txt"
#
# CnC = fromGI(C, fileC, colorPalette)
#
# ComponentsC = {
#     _:[]
#     for _ in range(nx.number_connected_components(C))
# }
# gIdx = 0
# tmpOldComp = []
# for node in C.nodes():
#     tmpNewComp = set(nx.node_connected_component(C, node))
#     if not tmpNewComp in tmpOldComp:
#         tmpOldComp.append(tmpNewComp)
#         ComponentsC[gIdx] = tmpNewComp
#         gIdx += 1
#
# print("You have an C level graph of %i nodes, %i edges, %i components." % (
#     len(C.nodes),
#     len(C.edges),
#     nx.number_connected_components(C)
# ))
#
#
#
# pickle.dump(ComponentsC,open(home + "/Data/outData/ComponentsCData.p", "wb" ))
# pickle.dump(ComponentsT,open(home + "/Data/outData/ComponentsTData.p", "wb" ))
# pickle.dump(ComponentsG,open(home + "/Data/outData/ComponentsGData.p", "wb" ))

# pickle.dump(C,open(home + "/Data/tmpData/GraphsCData.p", "wb" ))
# pickle.dump(T,open(home + "/Data/tmpData/GraphsTData.p", "wb" ))
# pickle.dump(G,open(home + "/Data/tmpData/GraphsGData.p", "wb" ))
