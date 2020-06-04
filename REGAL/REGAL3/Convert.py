


ethFile=open(f"{dataRoot}/REGAL/EtOH.edgelist",'wb')
nx.write_edgelist(ethG, ethFile, data=False, delimiter='\t')
ethFile.close()
dhtFile=open(f"{dataRoot}/REGAL/DHT.edgelist",'wb')
nx.write_edgelist(dhtG, dhtFile, data=False, delimiter='\t')
dhtFile.close()


nx.adjacency_matrix(dhtG).todense()
