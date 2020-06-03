

eth=/kuacc/users/ualtintas20/github/Data/CisGraph/Vers1.0/REGAL/EtOH.edgelist
dht=/kuacc/users/ualtintas20/github/Data/CisGraph/Vers1.0/REGAL/DHT.edgelist

bash CodeAndTestData/list2leda.sh $eth >> EtOH.gw
bash CodeAndTestData/list2leda.sh $dht >> DHT.gw

python ./GRAALRunner.py 0.8 DHT.gw EtOH.gw result
