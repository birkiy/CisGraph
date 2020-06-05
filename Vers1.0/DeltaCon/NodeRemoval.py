

from Function.Packages import *

ethG = pickle.load(open(f"{dataRoot}/tmpData/GraphsG.EtOH.Data.p", "rb" ))
dhtG = pickle.load(open(f"{dataRoot}/tmpData/GraphsG.DHT.Data.p", "rb" ))

dNodes = list(dhtG.nodes())


DF = pd.DataFrame()
for node in dNodes:
    df = pd.DataFrame()
    G = dhtG
    G.remove_node(node)
    df["node"] = node
    df["nodeClass"] = dhtG.nodes[node]["nodeClass"]
    df["Connected Components"] = nx.number_connected_components(G)
    DF = pd.concat([DF,df])
