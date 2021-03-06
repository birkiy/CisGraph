



groDmsMin=/home/ualtintas/github/Data/CisGraph/Vers2.0/Features/GroSeq/groseq.dmso.-.BPM.bigWig
groDmsPlu=/home/ualtintas/github/Data/CisGraph/Vers2.0/Features/GroSeq/groseq.dmso.+.BPM.bigWig
groDhtMin=/home/ualtintas/github/Data/CisGraph/Vers2.0/Features/GroSeq/groseq.dht.-.BPM.bigWig
groDhtPlu=/home/ualtintas/github/Data/CisGraph/Vers2.0/Features/GroSeq/groseq.dht.+.BPM.bigWig

heatmap=/home/ualtintas/github/Data/CisGraph/Vers2.0/Heatmaps
groPwd=/home/ualtintas/github/Data/CisGraph/Vers2.0/Features/GroSeq/groseq
homePwd=/home/ualtintas/github/Data/CisGraph/Vers1.0/Gro

# Promoters
tsP=~/genomeAnnotations/Regions/TSS.hg19.+.bed
tsM=~/genomeAnnotations/Regions/TSS.hg19.-.bed
tss=~/genomeAnnotations/Regions/TSS.hg19.Idx.bed
pro=~/ARBSs/regions/promoters_ann_5kb.bed

# ARBSs

oth=~/ARBSs/regions/otherARBS.bed
con=~/ARBSs/regions/cons-arbs.bed
ind=~/ARBSs/regions/ind-arbs.bed
non=~/ARBSs/regions/Non-Active-ARBS.bed
nAR=~/ARBSs/regions/negativeControl.ARBS.bed


echo "TSS centric CRE"

computeMatrix reference-point -S $groDmsMin $groDmsPlu $groDhtMin $groDhtPlu -R $oth $con $ind $non --skipZeros --referencePoint center -a 2000 -b 2000 -o $heatmap/ARBSMat.mat -p 30 -q

plotHeatmap -m $heatmap/ARBSMat.mat -out $heatmap/ARBSMatH.pdf --colorMap "YlGnBu" --whatToShow "heatmap and colorbar"

plotProfile --numPlotsPerRow 2 -m $heatmap/ARBSMat.mat -out $heatmap/ARBSMatP.pdf --color "#C0C0C0" "#5A5A5A" "#F9746D" "#ACACAC"




# FirstPart
awk -F'\t' '{print $1"\t"$2"\t"($3 - 350) " ""\t"$4
}' $tsP > tsP.FirstPart.bed

tsPF=~/genomeAnnotations/Regions/tsP.FirstPart.bed

awk -F'\t' '{print $1"\t"$2"\t"($3 - 350) " ""\t"$4
}' $tsM > tsM.FirstPart.bed

tsMF=~/genomeAnnotations/Regions/tsM.FirstPart.bed

awk -F'\t' '{print $1"\t"$2"\t"($3 - 350) " ""\t"$4
}' $tss > tss.FirstPart.bed

tssF=~/genomeAnnotations/Regions/tss.FirstPart.bed

awk -F'\t' '{print $1"\t"$2"\t"($3 - 350) " ""\t"$4
}' $con > con.FirstPart.bed

conF=~/ARBSs/regions/con.FirstPart.bed

awk -F'\t' '{print $1"\t"$2"\t"($3 - 350) " ""\t"$4
}' $ind > ind.FirstPart.bed

indF=~/ARBSs/regions/ind.FirstPart.bed

awk -F'\t' '{print $1"\t"$2"\t"($3 - 350) " ""\t"$4
}' $non > non.FirstPart.bed

nonF=~/ARBSs/regions/non.FirstPart.bed


awk -F'\t' '{print $1"\t"$2"\t"(int(($3 + $2)/2)) " ""\t"$4
}' $nAR > nAR.FirstPart.bed

nARF=~/ARBSs/regions/nAR.FirstPart.bed

awk -F'\t' '{print $1"\t"$2"\t"(int(($3 + $2)/2)) " ""\t"$4
}' $oth > oth.FirstPart.bed

othF=~/ARBSs/regions/oth.FirstPart.bed


# SecondPart
awk -F'\t' '{print $1"\t"($2 + 350) " ""\t"$3"\t"$4
}' $tsP > tsP.SecondPart.bed

tsPS=~/genomeAnnotations/Regions/tsP.SecondPart.bed
OvohF5el

awk -F'\t' '{print $1"\t"($2 + 350) " ""\t"$3"\t"$4
}' $tsM > tsM.SecondPart.bed

tsMS=~/genomeAnnotations/Regions/tsM.SecondPart.bed

awk -F'\t' '{print $1"\t"($2 + 350) " ""\t"$3"\t"$4
}' $tss > tss.SecondPart.bed

tssS=~/genomeAnnotations/Regions/tss.SecondPart.bed


awk -F'\t' '{print $1"\t"($2 + 350) " ""\t"$3"\t"$4
}' $con > con.SecondPart.bed

conS=~/ARBSs/regions/con.SecondPart.bed

awk -F'\t' '{print $1"\t"($2 + 350) " ""\t"$3"\t"$4
}' $ind > ind.SecondPart.bed

indS=~/ARBSs/regions/ind.SecondPart.bed

awk -F'\t' '{print $1"\t"($2 + 350) " ""\t"$3"\t"$4
}' $non > non.SecondPart.bed

nonS=~/ARBSs/regions/non.SecondPart.bed

awk -F'\t' '{print $1"\t"(int(($3 + $2)/2)) " ""\t"$3"\t"$4
}' $nAR > nAR.SecondPart.bed

nARS=~/ARBSs/regions/nAR.SecondPart.bed

awk -F'\t' '{print $1"\t"(int(($3 + $2)/2)) " ""\t"$3"\t"$4
}' $oth > oth.SecondPart.bed

othS=~/ARBSs/regions/oth.SecondPart.bed


cres=(
"con.FirstPart"
"con.SecondPart"
"ind.FirstPart"
"ind.SecondPart"
"non.FirstPart"
"non.SecondPart"
"nAR.FirstPart"
"nAR.SecondPart"
"oth.FirstPart"
"oth.SecondPart"
"tsP.FirstPart"
"tsP.SecondPart"
"tsM.FirstPart"
"tsM.SecondPart"
"tss.FirstPart"
"tss.SecondPart"
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
