

strTeth=/groups/lackgrp/ll_members/tunc/phd/ana-starrseq-lncap-lacklab/analysis/bigwig/lncap-logFC-etoh-plasmid.bigWig
strTdht=/groups/lackgrp/ll_members/tunc/phd/ana-starrseq-lncap-lacklab/analysis/bigwig/lncap-logFC-dht-plasmid.bigWig

creDhtData=/home/ualtintas/github/Data/CisGraph/Vers2.0/DHS/DHT_peaks_tabbed.bed
creEthData=/home/ualtintas/github/Data/CisGraph/Vers2.0/DHS/EtOH_peaks_tabbed.bed

feaData=/home/ualtintas/github/Data/CisGraph/Vers2.0/Features
heatmap=/home/ualtintas/github/Data/CisGraph/Vers2.0/Heatmaps

starrseq=$feaData/StarrSeq
groseq=$feaData/GroSeq
cage=$feaData/Cage

strEth=$starrseq/starrseq.etoh.plasmid.BPM.bigWig
strDht=$starrseq/starrseq.dht.plasmid.BPM.bigWig


groDmsMin=$groseq/groseq.dmso.BPM.minus.bigWig
groDmsPlu=$groseq/groseq.dmso.BPM.plus.bigWig
groDhtMin=$groseq/groseq.dht.BPM.minus.bigWig
groDhtPlu=$groseq/groseq.dht.BPM.plus.bigWig

tsP=/home/ualtintas/genomeAnnotations/Regions/TSS.hg19.+.bed
tsM=/home/ualtintas/genomeAnnotations/Regions/TSS.hg19.-.bed


con=/home/ualtintas/ARBSs/regions/cons-arbs.bed
ind=/home/ualtintas/ARBSs/regions/ind-arbs.bed
non=/home/ualtintas/ARBSs/regions/Non-Active-ARBS.bed


computeMatrix reference-point -S $strDht $strEth -R $creDhtData $creEthData --skipZeros --referencePoint center -a 2000 -b 2000 -o $heatmap/creStarrMat.mat -p 30 -q

plotHeatmap -m $heatmap/creStarrMat.mat -out $heatmap/cerStarrMat.pdf --heatmapWidth 10 --heatmapHeight 50 --colorMap "YlGnBu"

#

computeMatrix reference-point -S $strDht $strEth -R $creDhtData $creEthData --skipZeros --referencePoint center -a 2000 -b 2000 -o $heatmap/creStarrMatE.mat -p 30 -q

plotHeatmap -m $heatmap/creStarrMatE.mat -out $heatmap/cerStarrMatE.pdf --heatmapWidth 10 --heatmapHeight 50 --colorMap "YlGnBu"



computeMatrix reference-point -S $strDht $strEth -R $con $ind $non --skipZeros --referencePoint center -a 2000 -b 2000 -o $heatmap/arbsStarrMat.mat -p 30 -q

plotHeatmap -m $heatmap/arbsStarrMat.mat -out $heatmap/arbsStarrMat.pdf --heatmapWidth 10 --heatmapHeight 50 --colorMap "YlGnBu"




# computeMatrix reference-point -S $strDht $strEth -R $tsP $tsM --skipZeros --referencePoint center -a 2000 -b 2000 -o $heatmap/creStarrMat.mat -p 30 -q
#
# plotHeatmap -m $heatmap/creStarrMat.mat -out $heatmap/cerStarrMat.pdf --heatmapWidth 10 --heatmapHeight 50 --colorMap "YlGnBu"
#
#
#
# computeMatrix reference-point -S $strDht $strEth $strTeth $strTdht -R $creDhtData $creEthData --skipZeros --referencePoint center -a 2000 -b 2000 -o $heatmap/creStarrMat2.mat -p 30 -q
#
# plotHeatmap -m $heatmap/creStarrMat2.mat -out $heatmap/cerStarrMat2.pdf --heatmapWidth 10 --heatmapHeight 50 --colorMap "YlGnBu"





computeMatrix reference-point -S $groDmsMin $groDmsPlu $groDhtMin $groDhtPlu -R $cage/creCtr.bed $cage/creCtrP.bed $cage/creCtrM.bed --skipZeros --referencePoint center -a 2000 -b 2000 -o $heatmap/creGroCMat.mat -p 30 -q

plotHeatmap -m $heatmap/creGroCMat.mat -out $heatmap/creGroCMat.pdf --heatmapWidth 10 --heatmapHeight 50 --colorMap "YlGnBu"



computeMatrix reference-point -S $groDmsMin $groDmsPlu $groDhtMin $groDhtPlu -R $cage/creCtr.bed --skipZeros --referencePoint center -a 2000 -b 2000 -o $heatmap/creGroCBMat.mat -p 30 -q

plotHeatmap -m $heatmap/creGroCBMat.mat -out $heatmap/creGroCBMat.pdf --heatmapWidth 10 --heatmapHeight 50 --colorMap "YlGnBu"


computeMatrix reference-point -S $groDmsMin $groDmsPlu $groDhtMin $groDhtPlu -R $cage/creCtrP.bed $cage/creCtrM.bed --skipZeros --referencePoint center -a 2000 -b 2000 -o $heatmap/creGroCDMat.mat -p 30 -q

plotHeatmap -m $heatmap/creGroCDMat.mat -out $heatmap/creGroCDMat.pdf --heatmapWidth 10 --heatmapHeight 50 --colorMap "YlGnBu"




computeMatrix scale-regions -S $groDmsMin $groDmsPlu $groDhtMin $groDhtPlu -R $cage/creCtr.bed --skipZeros --regionBodyLength 500 -a 2000 -b 2000 -o $heatmap/creGroCSMat.mat -p 30 -q

plotHeatmap -m $heatmap/creGroCSMat.mat -out $heatmap/creGroCSMat.pdf --heatmapWidth 10 --heatmapHeight 50 --colorMap "YlGnBu"




computeMatrix reference-point -S $groDmsMin $groDmsPlu $groDhtMin $groDhtPlu -R $creDhtData $creEthData --skipZeros --referencePoint center -a 2000 -b 2000 -o $heatmap/creGroMat.mat -p 30 -q

plotHeatmap -m $heatmap/creGroMat.mat -out $heatmap/creGroMat.pdf --heatmapWidth 10 --heatmapHeight 50 --colorMap "YlGnBu"

####,





computeMatrix reference-point -S $groDmsMin $groDmsPlu $groDhtMin $groDhtPlu -R $tsP $tsM --skipZeros --referencePoint center -a 2000 -b 2000 -o $heatmap/tssGroMat.mat -p 30 -q

plotHeatmap -m $heatmap/tssGroMat.mat -out $heatmap/tssGroMat.pdf --heatmapWidth 10 --heatmapHeight 50 --colorMap "YlGnBu"



computeMatrix reference-point -S $cagEthPlu $cagEthMin $cagDhtPlu $cagDhtMin -R $creDhtData $creEthData --skipZeros --referencePoint center -a 2000 -b 2000 -o $heatmap/creCageMat.mat -p 30 -q

plotHeatmap -m $heatmap/creCageMat.mat -out $heatmap/creCageMat.pdf --heatmapWidth 10 --heatmapHeight 50 --colorMap "YlGnBu"




computeMatrix reference-point -S $cagEthPlu $cagEthMin $cagDhtPlu $cagDhtMin -R $tsP $tsM --skipZeros --referencePoint center -a 2000 -b 2000 -o $heatmap/tssCageMat.mat -p 30 -q

plotHeatmap -m $heatmap/tssCageMat.mat -out $heatmap/tssCageMat.pdf --heatmapWidth 10 --heatmapHeight 50 --colorMap "YlGnBu"




computeMatrix reference-point -S $cage/*.bw -R $tsP $tsM --skipZeros --referencePoint center -a 2000 -b 2000 -o $heatmap/tssCageMat.mat -p 30 -q

plotHeatmap -m $heatmap/tssCageMat.mat -out $heatmap/tssCageMat.pdf --heatmapWidth 10 --heatmapHeight 50 --colorMap "YlGnBu"





computeMatrix reference-point -S \
$outData/LNCaP.0h.rep1.-.BPM.bw \
$outData/LNCaP.0h.rep1.+.BPM.bw \
$outData/LNCaP.3h.rep1.-.BPM.bw \
$outData/LNCaP.3h.rep1.+.BPM.bw \
$outData/LNCaP.6h.rep1.-.BPM.bw \
$outData/LNCaP.6h.rep1.+.BPM.bw \
$outData/LNCaP.18h.rep1.-.BPM.bw \
$outData/LNCaP.18h.rep1.+.BPM.bw \
$outData/LNCaP.48h.rep1.-.BPM.bw \
$outData/LNCaP.48h.rep1.+.BPM.bw \
$outData/LNCaP.72h.rep1.-.BPM.bw \
$outData/LNCaP.72h.rep1.+.BPM.bw \
-R $tsP $tsM --skipZeros --referencePoint center -a 2000 -b 2000 -o $heatmap/tssCageMat.mat -p 30 -q

plotHeatmap -m $heatmap/tssCageMat.mat -out $heatmap/tssCageMat.pdf --heatmapWidth 10 --heatmapHeight 50 --colorMap "YlGnBu"

computeMatrix reference-point -S \
$outData/LNCaP.0h.rep1.-.BPM.bw \
$outData/LNCaP.0h.rep1.+.BPM.bw \
$outData/LNCaP.3h.rep1.-.BPM.bw \
$outData/LNCaP.3h.rep1.+.BPM.bw \
$outData/LNCaP.6h.rep1.-.BPM.bw \
$outData/LNCaP.6h.rep1.+.BPM.bw \
$outData/LNCaP.18h.rep1.-.BPM.bw \
$outData/LNCaP.18h.rep1.+.BPM.bw \
$outData/LNCaP.48h.rep1.-.BPM.bw \
$outData/LNCaP.48h.rep1.+.BPM.bw \
$outData/LNCaP.72h.rep1.-.BPM.bw \
$outData/LNCaP.72h.rep1.+.BPM.bw \
-R $creDhtData $creEthData --skipZeros --referencePoint center -a 500 -b 500 -o $heatmap/creCageMat5.mat -p 30 -q --outFileNameMatrix creCageMat5.tab --outFileSortedRegions creCageMat5.bed


plotHeatmap -m $heatmap/creCageMat5.mat -out $heatmap/creCageMat5.pdf --heatmapWidth 10 --heatmapHeight 50 --colorMap "YlGnBu"


tail -n +4 creCageMat5.tab > headles.tab
sed -n '3p' creCageMat5.tab | cut -f 3- > header
cat header headles.tab > creCageMat.modified.tab
paste creCageMat5.bed creCageMat.modified.tab > creCageConcat.tab







computeMatrix reference-point -S \
$outData/LNCaP.0h.rep1.-.BPM.bw \
$outData/LNCaP.0h.rep1.+.BPM.bw \
$outData/LNCaP.3h.rep1.-.BPM.bw \
$outData/LNCaP.3h.rep1.+.BPM.bw \
$outData/LNCaP.6h.rep1.-.BPM.bw \
$outData/LNCaP.6h.rep1.+.BPM.bw \
$outData/LNCaP.18h.rep1.-.BPM.bw \
$outData/LNCaP.18h.rep1.+.BPM.bw \
$outData/LNCaP.48h.rep1.-.BPM.bw \
$outData/LNCaP.48h.rep1.+.BPM.bw \
$outData/LNCaP.72h.rep1.-.BPM.bw \
$outData/LNCaP.72h.rep1.+.BPM.bw \
-R $con $ind $non --skipZeros --referencePoint center -a 350 -b 350 -o $heatmap/ArbsCageMat35.mat -p 30 -q --outFileNameMatrix ArbsCageMat35.tab --outFileSortedRegions ArbsCageMat35.bed


plotHeatmap -m $heatmap/ArbsCageMat35.mat -out $heatmap/ArbsCageMat35.pdf --heatmapWidth 10 --heatmapHeight 50 --colorMap "YlGnBu"




computeMatrix reference-point -S \
$outData/cageMerged.-.BPM.bw \
$outData/cageMerged.+.BPM.bw \
-R $con $ind $non --skipZeros --referencePoint center -a 350 -b 350 -o $heatmap/ArbsCageMMat35.mat -p 30 -q

plotHeatmap -m $heatmap/ArbsCageMMat35.mat -out $heatmap/ArbsCageMMat35.pdf --heatmapWidth 10 --heatmapHeight 50 --colorMap "YlGnBu"





computeMatrix reference-point -S $cage/*.bigWig -R $tsP $tsM --skipZeros --referencePoint center -a 2000 -b 2000 -o $heatmap/tssCageB2BMat.mat -p 30 -q

plotHeatmap -m $heatmap/tssCageB2BMat.mat -out $heatmap/tssCageB2BMat.pdf --heatmapWidth 10 --heatmapHeight 50 --colorMap "YlGnBu"



computeMatrix reference-point -S $cage/*raw* -R $tsP $tsM --skipZeros --referencePoint center -a 2000 -b 2000 -o $heatmap/tssCageRMat.mat -p 30 -q

plotHeatmap -m $heatmap/tssCageRMat.mat -out $heatmap/tssCageRMat.pdf --heatmapWidth 10 --heatmapHeight 50 --colorMap "YlGnBu"



computeMatrix reference-point -S $cage/*24h* $cage/*0h*  -R $creDhtData $creEthData --skipZeros --referencePoint center -a 2000 -b 2000 -o $heatmap/creCageRMat.mat -p 30 -q

plotHeatmap -m $heatmap/creCageRMat.mat -out $heatmap/creCageRMat.pdf --heatmapWidth 10 --heatmapHeight 50 --colorMap "YlGnBu"




computeMatrix reference-point -S $cage/*raw* -R $con $ind $non --skipZeros --referencePoint center -a 2000 -b 2000 -o $heatmap/ArbsCageRMat.mat -p 30 -q

plotHeatmap -m $heatmap/ArbsCageRMat.mat -out $heatmap/ArbsCageRMat.pdf --heatmapWidth 10 --heatmapHeight 50 --colorMap "YlGnBu"




computeMatrix reference-point -S $cage/*raw* -R $tsP $tsM --skipZeros --referencePoint center -a 2000 -b 2000 -o $heatmap/tssCageRMat.mat -p 30 -q

plotHeatmap -m $heatmap/tssCageRMat.mat -out $heatmap/tssCageRMat.pdf --heatmapWidth 10 --heatmapHeight 50 --colorMap "YlGnBu"








computeMatrix reference-point -S $cagEthPlu $cagEthMin $cagDhtPlu $cagDhtMin -R $con $ind $non --skipZeros --referencePoint center -a 10000 -b 10000 -o $heatmap/arbsCageMat.mat -p 30 -q

plotHeatmap -m $heatmap/arbsCageMat.mat -out $heatmap/arbsCageMat.pdf --heatmapWidth 10 --heatmapHeight 50 --colorMap "YlGnBu"





computeMatrix reference-point -S $cagEthPlu $cagEthMin $cagDhtPlu $cagDhtMin $strDht $strEth $groDmsMin $groDmsPlu $groDhtMin $groDhtPlu -R $creDhtData $creEthData --skipZeros --referencePoint center -a 2000 -b 2000 -o $heatmap/creMat.mat -p 30 -q

plotHeatmap -m $heatmap/creMat.mat -out $heatmap/cerMat.pdf --heatmapWidth 10 --heatmapHeight 50 --colorMap "YlGnBu"
