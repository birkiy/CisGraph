

dht=~/github/Data/CisGraph/Vers1.0/InteractionBedPe/LNCaP_DHT_2000_CiceroConns.bedpe

eth=~/github/Data/CisGraph/Vers1.0/InteractionBedPe/LNCaP_EtOH_2000_CiceroConns.bedpe


awk -F'\t' '{if($7 != "NA"){print}}' $dht > LNCaP_DHT_2000_CiceroConns.woNA.bedpe


cat LNCaP_DHT_2000_CiceroConns.woNA.bedpe |tr ',' '.' > LNCaP_DHT_2000_CiceroConns.woNA.dot.bedpe


dht=LNCaP_DHT_2000_CiceroConns.woNA.dot.bedpe


awk -F'\t' '{if($7 != "NA"){print}}' $eth > LNCaP_EtOH_2000_CiceroConns.woNA.bedpe


cat LNCaP_EtOH_2000_CiceroConns.woNA.bedpe |tr ',' '.' > LNCaP_EtOH_2000_CiceroConns.woNA.dot.bedpe


eth=LNCaP_EtOH_2000_CiceroConns.woNA.dot.bedpe

# Ultimately, Cicero provides a "Cicero co-accessibility" score between -1 and 1
# between each pair of accessibile peaks within a user defined distance where a
# higher number indicates higher co-accessibility.
awk -F'\t' '($7 + 0) > 0.25' $dht > LNCaP_DHT_2000_CiceroConns.25Filter.bedpe
awk -F'\t' '($7 + 0) > 0.20' $dht > LNCaP_DHT_2000_CiceroConns.20Filter.bedpe
awk -F'\t' '($7 + 0) > 0.10' $dht > LNCaP_DHT_2000_CiceroConns.10Filter.bedpe
awk -F'\t' '($7 + 0) > 0.05' $dht > LNCaP_DHT_2000_CiceroConns.05Filter.bedpe




awk -F'\t' '($7 + 0) > 0.25' $eth > LNCaP_EtOH_2000_CiceroConns.25Filter.bedpe
awk -F'\t' '($7 + 0) > 0.20' $eth > LNCaP_EtOH_2000_CiceroConns.20Filter.bedpe
awk -F'\t' '($7 + 0) > 0.10' $eth > LNCaP_EtOH_2000_CiceroConns.10Filter.bedpe
awk -F'\t' '($7 + 0) > 0.05' $eth > LNCaP_EtOH_2000_CiceroConns.05Filter.bedpe






awk -F'\t' '{if($7 > 0.05){print}}' $dht > LNCaP_DHT_2000_CiceroConns.05Filter.bedpe

awk -F'\t' '{if($7 > 0.05){print}}' $eth > LNCaP_EtOH_2000_CiceroConns.05Filter.bedpe



awk -F'\t' '{if($7 > 0.01){print}}' $dht > LNCaP_DHT_2000_CiceroConns.01Filter.bedpe

awk -F'\t' '{if($7 > 0.01){print}}' $eth > LNCaP_EtOH_2000_CiceroConns.01Filter.bedpe



awk -F'\t' '{if($7 > 0.00){print}}' $dht > LNCaP_DHT_2000_CiceroConns.00Filter.bedpe

awk -F'\t' '{if($7 > 0.00){print}}' $eth > LNCaP_EtOH_2000_CiceroConns.00Filter.bedpe




awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {if(abs($7) > 0.1){print}}' $dht > LNCaP_DHT_2000_CiceroConns.10Filter.bedpe

awk -F'\t' '{if($7 > 0.1){print}}' $eth > LNCaP_EtOH_2000_CiceroConns.10Filter.bedpe




awk -F'\t' '{if($7 > 0.0001){print}}' $dht > LNCaP_DHT_2000_CiceroConns.0001Filter.bedpe

awk -F'\t' '{if($7 > 0.0001){print}}' $eth > LNCaP_EtOH_2000_CiceroConns.0001Filter.bedpe



awk -F'\t' '{if(($7) > 0.25){print}}' $dht > LNCaP_DHT_2000_CiceroConns.25Filter.bedpe

awk -F'\t' '{if(($7) > 0.25){print}}' $eth > LNCaP_EtOH_2000_CiceroConns.25Filter.bedpe






awk -F'\t' '{if($7 > 0.00001){print}}' LNCaP_DHT_2000_CiceroConns.bedpe > LNCaP_DHT_2000_CiceroConns.01Filter.bedpe

awk -F'\t' '{if($7 > 0.00001){print}}' $eth > LNCaP_EtOH_2000_CiceroConns.00001Filter.bedpe








647436

896990
