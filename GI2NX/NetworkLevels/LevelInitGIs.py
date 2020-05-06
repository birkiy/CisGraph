

from Functions.GIFunctions import *





colorPalette = {"enh": "#a3c9a8",

                "con": "#FADE89",
                "ind": "#57A4B1",
                "non": "#B0D894",

                "arp": "#F99FA5",
                "ctc": "#9FE2D4",
                "fox": "#C1A2EF",
                "med": "#9ECCF2",
                "enz": "#F7C29C",
                "ari": "#F9E0A1",
                "pol": "#F29ECC",
                "nmy": "#E2B79D",

                "gen": "#ddd8c4",

                "upP": "#95342c",
                "dwP": "#2c6589",

                "tad": "#84b59f",
                "com": "#69a297" ,
                "chr": "#50808e"
}





G = nx.Graph()
fileG = home + "/Data/GIs/GI.G.txt"

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
fileT = home + "/Data/GIs/GI.T.txt"

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
fileC = home + "/Data/GIs/GI.C.txt"

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



pickle.dump(ComponentsC,open(home + "/Data/outData/ComponentsCData.p", "wb" ))
pickle.dump(ComponentsT,open(home + "/Data/outData/ComponentsTData.p", "wb" ))
pickle.dump(ComponentsG,open(home + "/Data/outData/ComponentsGData.p", "wb" ))

pickle.dump(C,open(home + "/Data/tmpData/GraphsCData.p", "wb" ))
pickle.dump(T,open(home + "/Data/tmpData/GraphsTData.p", "wb" ))
pickle.dump(G,open(home + "/Data/tmpData/GraphsGData.p", "wb" ))
