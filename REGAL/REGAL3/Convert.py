


ethFile=open(f"{dataRoot}/REGAL/EtOH.edgelist",'wb')
nx.write_edgelist(ethG, ethFile, data=False)

dhtFile=open(f"{dataRoot}/REGAL/DHT.edgelist",'wb')
nx.write_edgelist(dhtG, dhtFile, data=False)

nx.adjacency_matrix(dhtG).todense()
