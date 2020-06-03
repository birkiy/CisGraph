

dht=~/github/Data/CisGraph/Vers1.0/InteractionBedPe/LNCaP_DHT_2000_CiceroConns.bedpe

eth=~/github/Data/CisGraph/Vers1.0/InteractionBedPe/LNCaP_EtOH_2000_CiceroConns.bedpe

awk -F'\t' '{if($7 > 0.25){print}}' $dht > LNCaP_DHT_2000_CiceroConns.25Filter.bedpe

awk -F'\t' '{if($7 > 0.25){print}}' $eth > LNCaP_EtOH_2000_CiceroConns.25Filter.bedpe
