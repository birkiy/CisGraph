
# password = OvohF5el


min=/groups/lackgrp/ll_members/common/LNCaP_bigwigs/groseq_analysis/mapping/groSeq_lncap_minus_dht-dmso.bigWig
plu=/groups/lackgrp/ll_members/common/LNCaP_bigwigs/groseq_analysis/mapping/groSeq_lncap_plus_dht-dmso.bigWig

con=~/regions/cons-arbs.bed
ind=~/regions/ind-arbs.bed
non=~/regions/Non-Active-ARBS.bed
pro=~/regions/promoters_ann_5kb.bed

uPP=~/genomeAnnotations/TSS.hg19.+.up.bed
dPP=~/genomeAnnotations/TSS.hg19.+.dw.bed

uMP=~/genomeAnnotations/TSS.hg19.-.up.bed
dMP=~/genomeAnnotations/TSS.hg19.-.dw.bed


minDht=/groups/lackgrp/ll_members/common/LNCaP_bigwigs/groseq_analysis/mapping/groseq_dht.minus.bw
pluDht=/groups/lackgrp/ll_members/common/LNCaP_bigwigs/groseq_analysis/mapping/groseq_dht.plus.bw


minDms=/groups/lackgrp/ll_members/common/LNCaP_bigwigs/groseq_analysis/mapping/groseq_dmso.minus.bw
pluDms=/groups/lackgrp/ll_members/common/LNCaP_bigwigs/groseq_analysis/mapping/groseq_dmso.plus.bw


computeMatrix reference-point -S $min $plu -R $con $ind $non $pro -a 3000 -b 3000 -o groRegions.mat.gz -p 10
computeMatrix reference-point -S $min $plu -R $con $ind $non -a 3000 -b 3000 -o groARBSRegions.mat.gz -p 10


computeMatrix reference-point -S $min $plu -R $con $ind $non -a 20000 -b 20000 -o groARBSRegions+-20.mat.gz -p 10

plotHeatmap -m groRegions.mat.gz -out groHeatmap.pdf --heatmapWidth 10 --heatmapHeight 70 --colorMap "YlGnBu"

plotHeatmap -m groARBSRegions.mat.gz -out groARBSHeatmap.pdf --heatmapWidth 10 --heatmapHeight 70 --colorMap "YlGnBu"

plotHeatmap -m groARBSRegions+-20.mat.gz -out groARBS+-20Heatmap.pdf --heatmapWidth 10 --heatmapHeight 70 --colorMap "YlGnBu"







computeMatrix reference-point -S $min $plu -R $con -a 3000 -b 3000 -o groCon.mat.gz -p 15

computeMatrix reference-point -S $min $plu -R $ind -a 3000 -b 3000 -o groInd.mat.gz -p 15

computeMatrix reference-point -S $min $plu -R $non -a 3000 -b 3000 -o groNot.mat.gz -p 15


plotHeatmap -m groCon.mat.gz -out conHeatmap.pdf --heatmapWidth 10 --heatmapHeight 70 --colorMap "YlGnBu" --hclust 2

plotHeatmap -m groInd.mat.gz -out indHeatmap.pdf --heatmapWidth 10 --heatmapHeight 70 --colorMap "YlGnBu" --hclust 2

plotHeatmap -m groNon.mat.gz -out nonHeatmap.pdf --heatmapWidth 10 --heatmapHeight 70 --colorMap "YlGnBu" --hclust 2






computeMatrix reference-point -S $minDht $pluDht -R $con $ind $non -a 3000 -b 3000 -o groARBSRegionsDHT.mat.gz -p 15

computeMatrix reference-point -S $minDms $pluDms -R $con $ind $non -a 3000 -b 3000 -o groARBSRegionsDMSO.mat.gz -p 15

plotHeatmap -m groARBSRegionsDHT.mat.gz -out groHeatmapDHT.pdf --heatmapWidth 10 --heatmapHeight 70 --colorMap "YlGnBu"

plotHeatmap -m groARBSRegionsDMSO.mat.gz -out groHeatmapDMSO.pdf --heatmapWidth 10 --heatmapHeight 70 --colorMap "YlGnBu"








computeMatrix reference-point -S $minDht $pluDht $minDms $pluDms -R $upP $dwP $con $ind $non -a 3000 -b 3000 -o groARBSProRegionsDHTvsDMSO.mat.gz -p 15

plotHeatmap -m groARBSProRegionsDHTvsDMSO.mat.gz -out groHeatmapProARBSDHTvsDMSO.pdf --heatmapWidth 10 --heatmapHeight 70 --colorMap "YlGnBu"




computeMatrix reference-point -S $minDht $pluDht $minDms $pluDms -R $con $ind $non -a 3000 -b 3000 -o groARBSRegionsDHTvsDMSO.mat.gz -p 15

plotHeatmap -m groARBSRegionsDHTvsDMSO.mat.gz -out groHeatmapARBSDHTvsDMSO.pdf --heatmapWidth 10 --heatmapHeight 70 --colorMap "YlGnBu"





computeMatrix reference-point -S $minDht $pluDht $minDms $pluDms -R $pro -a 3000 -b 3000 -o groProRegionsDHTvsDMSO.mat.gz -p 15

plotHeatmap -m groProRegionsDHTvsDMSO.mat.gz -out groHeatmapProDHTvsDMSO.pdf --heatmapWidth 10 --heatmapHeight 70 --colorMap "YlGnBu" --hclust 5

# HOMER annotatePeaks
# # http://homer.ucsd.edu/homer/ngs/annotation.html
# TSS (by default defined from -1kb to +100bp)
# TTS (by default defined from -100 bp to +1kb)
annotatePeaks.pl tss hg19 -size -500,500 -cTSS > TSS.hg19.gtf
annotatePeaks.pl tts hg19 -size -500,500 -cTSS > TTS.hg19.gtf

cut TSS.hg19.gtf -f2,3,4,16 > tss.pos
# cut TSS.hg19.gtf -f1 > tss.id
cut TSS.hg19.gtf -f5 > tss.st

paste tss.pos tss.st > TSS.hg19.bed

awk -F'\t' '{if($5 == "+") {print}}' TSS.hg19.bed > TSS.hg19.+.bed
awk -F'\t' '{if($5 == "-") {print}}' TSS.hg19.bed > TSS.hg19.-.bed

# bedtools sort -i TSS.hg19.+.bed > TSS.hg19.+Sorted.bed
# bedtools sort -i TSS.hg19.-.bed > TSS.hg19.-Sorted.bed


cut TTS.hg19.gtf -f2,3,4,16 > tts.pos
# cut TTS.hg19.gtf -f1 > tts.id
cut TTS.hg19.gtf -f5 > tts.st

paste tts.pos tts.st > TTS.hg19.bed

awk -F'\t' '{if($5 == "+") {print}}' TTS.hg19.bed > TTS.hg19.+.bed
awk -F'\t' '{if($5 == "-") {print}}' TTS.hg19.bed > TTS.hg19.-.bed

# bedtools sort -i TTS.hg19.+.bed > TTS.hg19.+Sorted.bed
# bedtools sort -i TTS.hg19.-.bed > TTS.hg19.-Sorted.bed

# prP=~/genomeAnnotations/TSS.hg19.+.bed
# prM=~/genomeAnnotations/TSS.hg19.-.bed
#
# tsP=~/genomeAnnotations/TTS.hg19.+.bed
# tsM=~/genomeAnnotations/TTS.hg19.-.bed

cat ~/genomeAnnotations/TSS.hg19.+.dw.bed ~/genomeAnnotations/TSS.hg19.+.up.bed > TSS.hg19.+.DE.bed
cat ~/genomeAnnotations/TSS.hg19.-.dw.bed ~/genomeAnnotations/TSS.hg19.-.up.bed > TSS.hg19.-.DE.bed

prP=~/genomeAnnotations/TSS.hg19.+.DE.bed
prM=~/genomeAnnotations/TSS.hg19.-.DE.bed


computeMatrix reference-point -q -S $pluDht $minDht $pluDms $minDms -R $prP $prM $con $ind $non --referencePoint center -a 3000 -b 3000 -o groARBSProRegionsDHTvsDMSO.mat.gz -p 30

plotHeatmap -m groARBSProRegionsDHTvsDMSO.mat.gz -out groHeatmapARBSProDHTvsDMSO.pdf --heatmapWidth 10 --heatmapHeight 70 --colorMap "YlGnBu"


computeMatrix reference-point -q -S $pluDht $minDht $pluDms $minDms -R $con $ind $non --referencePoint center -a 3000 -b 3000 -o groARBSRegionsDHTvsDMSO.mat.gz -p 30

plotHeatmap -m groARBSRegionsDHTvsDMSO.mat.gz -out groHeatmapARBSDHTvsDMSO.pdf --heatmapWidth 10 --heatmapHeight 70 --colorMap "YlGnBu"





computeMatrix scale-regions -S $pluDht $minDht $pluDms $minDms -R $prP $prM $tsP $tsM -b 3000 -o groProRegionsDHTvsDMSO.scale.mat.gz -p 15

plotHeatmap -m groProRegionsDHTvsDMSO.scale.mat.gz -out groHeatmapProDHTvsDMSO.scale.pdf --heatmapWidth 10 --heatmapHeight 70 --colorMap "YlGnBu"





#################3

Homer TSS and TTS
