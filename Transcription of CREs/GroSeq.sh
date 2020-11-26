


groDmsMin=/home/ualtintas/github/Data/CisGraph/Vers2.0/Features/GroSeq/groseq.dmso.-.BPM.bigWig
groDmsPlu=/home/ualtintas/github/Data/CisGraph/Vers2.0/Features/GroSeq/groseq.dmso.+.BPM.bigWig
groDhtMin=/home/ualtintas/github/Data/CisGraph/Vers2.0/Features/GroSeq/groseq.dht.-.BPM.bigWig
groDhtPlu=/home/ualtintas/github/Data/CisGraph/Vers2.0/Features/GroSeq/groseq.dht.+.BPM.bigWig

groPwd=/home/ualtintas/github/Data/CisGraph/Vers2.0/Features/GroSeq/groseq

tss=~/genomeAnnotations/Regions/TSS.hg19.Idx.bed

con=~/ARBSs/regions/cons-arbs.bed
ind=~/ARBSs/regions/ind-arbs.bed
non=~/ARBSs/regions/Non-Active-ARBS.bed
nAR=~/ARBSs/regions/negativeControl.ARBS.bed

outPwd=/home/ualtintas/github/Data/CisGraph/Vers1.0/Gro


homePwd=~/genomeAnnotations/Regions
cres=(
"TSS.hg19.Idx.woS"
)

cut -f 1,2,3,4 $homePwd/TSS.hg19.Idx.bed > $homePwd/TSS.hg19.Idx.woS.bed

homePwd=~/ARBSs/regions
cres=(
"cons-arbs"
"ind-arbs"
"Non-Active-ARBS"
"negativeControl.ARBS"
)

gros=(
".dmso.-"
".dmso.+"
".dht.-"
".dht.+"
)

for i in ${cres[@]}
do
  cut -f4 $homePwd/$i".bed" > "$i.tmp"
  for j in ${gros[@]}
  do
  bigWigAverageOverBed $groPwd$j".BPM.bigWig" $homePwd/$i".bed" $outPwd/$i$j".tab";
  sort -k 1 $outPwd/$i$j".tab" | cut -f4 | paste -d"\t" "$i.tmp" - > "$i.tab"
  cat "$i.tab" > "$i.tmp"
  rm $outPwd/$i$j".tab"
done
rm -f "$i.tmp"
done
