



groDmsMin=/home/ualtintas/github/Data/CisGraph/Vers2.0/Features/GroSeq/groseq.dmso.-.BPM.bigWig
groDmsPlu=/home/ualtintas/github/Data/CisGraph/Vers2.0/Features/GroSeq/groseq.dmso.+.BPM.bigWig
groDhtMin=/home/ualtintas/github/Data/CisGraph/Vers2.0/Features/GroSeq/groseq.dht.-.BPM.bigWig
groDhtPlu=/home/ualtintas/github/Data/CisGraph/Vers2.0/Features/GroSeq/groseq.dht.+.BPM.bigWig

groPwd=/home/ualtintas/github/Data/CisGraph/Vers2.0/Features/GroSeq/groseq
homePwd=/home/ualtintas/github/Data/CisGraph/Vers1.0/Gro

# Promoters
tsP=~/genomeAnnotations/Regions/TSS.hg19.+.bed
tsM=~/genomeAnnotations/Regions/TSS.hg19.-.bed
pro=~/ARBSs/regions/promoters_ann_5kb.bed

# ARBSs

con=~/ARBSs/regions/cons-arbs.bed
ind=~/ARBSs/regions/ind-arbs.bed
non=~/ARBSs/regions/Non-Active-ARBS.bed


# FirstPart
awk -F'\t' '{print $1"\t"$2"\t"($3 - 500) " ""\t"$4"\t"$5
}' $tsP > tsP.FirstPart.bed

tsPF=~/genomeAnnotations/Regions/tsP.FirstPart.bed

awk -F'\t' '{print $1"\t"$2"\t"($3 - 500) " ""\t"$4"\t"$5
}' $tsM > tsM.FirstPart.bed

tsMF=~/genomeAnnotations/Regions/tsM.FirstPart.bed

awk -F'\t' '{print $1"\t"$2"\t"($3 - 350) " ""\t"$4"\t"$5
}' $con > con.FirstPart.bed

conF=~/genomeAnnotations/Regions/con.FirstPart.bed

awk -F'\t' '{print $1"\t"$2"\t"($3 - 350) " ""\t"$4"\t"$5
}' $ind > ind.FirstPart.bed

indF=~/genomeAnnotations/Regions/ind.FirstPart.bed

awk -F'\t' '{print $1"\t"$2"\t"($3 - 350) " ""\t"$4"\t"$5
}' $non > non.FirstPart.bed

nonF=~/genomeAnnotations/Regions/non.FirstPart.bed


# SecondPart
awk -F'\t' '{print $1"\t"($2 + 500) " ""\t"$3"\t"$4"\t"$5
}' $tsP > tsP.SecondPart.bed

tsPS=~/genomeAnnotations/Regions/tsP.SecondPart.bed
OvohF5el

awk -F'\t' '{print $1"\t"($2 + 500) " ""\t"$3"\t"$4"\t"$5
}' $tsM > tsM.SecondPart.bed

tsMS=~/genomeAnnotations/Regions/tsM.SecondPart.bed

awk -F'\t' '{print $1"\t"($2 + 350) " ""\t"$3"\t"$4"\t"$5
}' $con > con.SecondPart.bed

conS=~/genomeAnnotations/Regions/con.SecondPart.bed

awk -F'\t' '{print $1"\t"($2 + 350) " ""\t"$3"\t"$4"\t"$5
}' $ind > ind.SecondPart.bed

indS=~/genomeAnnotations/Regions/ind.SecondPart.bed

awk -F'\t' '{print $1"\t"($2 + 350) " ""\t"$3"\t"$4"\t"$5
}' $non > non.SecondPart.bed

nonS=~/genomeAnnotations/Regions/non.SecondPart.bed




cres=(
"tsP.FirstPart"
"tsP.SecondPart"
"tsM.FirstPart"
"tsM.SecondPart"
"con.FirstPart"
"con.SecondPart"
"ind.FirstPart"
"ind.SecondPart"
"non.FirstPart"
"non.SecondPart"
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
  bigWigAverageOverBed $groPwd$j".BPM.bigWig" $homePwd/$i".bed" $homePwd/$i$j".tab";
  sort -k 1 $homePwd/$i$j".tab" | cut -f4 | paste -d"\t" "$i.tmp" - > "$i.tab"
  cat "$i.tab" > "$i.tmp"
  rm $homePwd/$i$j".tab"
done
rm -f "$i.tmp"
done
