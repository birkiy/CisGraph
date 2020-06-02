


from BFS.CiceroBFS import *


setDhtRoot = set(dhtRooted.keys())
setEthRoot = set(ethRooted.keys())

both = setDhtRoot.intersection(setEthRoot)

PopUpD = {"con": [], "ind": [], "non": []}
PopUpE = {"con": [], "ind": [], "non": []}
for root in both:
    dht = set(dhtRooted[root])
    eth = set(ethRooted[root])
    dhtU = dht.difference(eth)
    ethU = eth.difference(dht)
    for n in ethU:
        PopUpE[Gs[1].nodes[n]["nodeClass"]] += [1/len(eth)]
    for n in dhtU:
        PopUpD[Gs[0].nodes[n]["nodeClass"]] += [1/len(eth)]









ethG = pickle.load(open(f"{dataRoot}/tmpData/GraphsG.EtOH.Data.p", "rb" ))
dhtG = pickle.load(open(f"{dataRoot}/tmpData/GraphsG.DHT.Data.p", "rb" ))

Gs = [dhtG, ethG]


setDht = set(Gs[0].keys())
setEth = set(Gs[1].keys())
