

ARBS=~/github/Data/CisGraph/Vers1.0/Regions/GSE83860_GSM2219854_ChIPseq_LNCaP_AR_DHT_sorted.bed

con=~/github/Data/CisGraph/Vers1.0/Regions/cons-arbs.bed
ind=~/github/Data/CisGraph/Vers1.0/Regions/ind-arbs.bed
non=~/github/Data/CisGraph/Vers1.0/Regions/Non-Active-ARBS.bed



bedtools intersect -a $ARBS -b $con $ind $non -wb | cut -f11 | sort | uniq > ARBS.txt


cat $con $ind $non | cut -f4 | sort | uniq > Enh.txt

comm -3 Enh.txt ARBS.txt

grep "peakno_2366-ARBS" $con
grep "peakno_2366-ARBS" $ind # > chr3	15714862	15715562	peakno_2366-ARBS
grep "peakno_2366-ARBS" $non


bedtools intersect -a $ARBS -b $con $ind $non -v > otherARBS.bed
