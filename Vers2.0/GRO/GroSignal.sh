


groDmsMin=/home/ualtintas/github/Data/CisGraph/Vers2.0/Features/GroSeq/groseq.dmso.-.BPM.bigWig
groDmsPlu=/home/ualtintas/github/Data/CisGraph/Vers2.0/Features/GroSeq/groseq.dmso.+.BPM.bigWig
groDhtMin=/home/ualtintas/github/Data/CisGraph/Vers2.0/Features/GroSeq/groseq.dht.-.BPM.bigWig
groDhtPlu=/home/ualtintas/github/Data/CisGraph/Vers2.0/Features/GroSeq/groseq.dht.+.BPM.bigWig

groPwd=/home/ualtintas/github/Data/CisGraph/Vers2.0/Features/GroSeq/groseq
cagePwd=/home/ualtintas/github/Data/CisGraph/Vers2.0/Features/Cage
cres=(
"creCtrA1"
"creCtrA2"
"creCtrB1"
"creCtrB2"
"creCtrM1"
"creCtrM2"
"creCtrP1"
"creCtrP2"
)

gros=(
".dmso.-"
".dmso.+"
".dht.-"
".dht.+"
)

for i in ${cres[@]}
do
  cut -f4 $cagePwd/$i".bed" > "$i.tmp"
  for j in ${gros[@]}
  do
  bigWigAverageOverBed $groPwd$j".BPM.bigWig" $cagePwd/$i".bed" $cagePwd/$i$j".tab";
  sort -k 1 $cagePwd/$i$j".tab" | cut -f4 | paste -d"\t" "$i.tmp" - > "$i.tab"
  cat "$i.tab" > "$i.tmp"
  rm $cagePwd/$i$j".tab"
done
rm -f "$i.tmp"
done


python GroDirectionSignal.py
