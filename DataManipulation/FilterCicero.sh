

dht=~/github/Data/CisGraph/Vers1.0/InteractionBedPe/LNCaP_DHT_2000_CiceroConns.bedpe

eth=~/github/Data/CisGraph/Vers1.0/InteractionBedPe/LNCaP_EtOH_2000_CiceroConns.bedpe

awk -F'\t' '{if($7 > 0.05){print}}' $dht > LNCaP_DHT_2000_CiceroConns.Modified.bedpe

awk -F'\t' '{if($7 > 0.05){print}}' $eth > LNCaP_EtOH_2000_CiceroConns.Modified.bedpe
