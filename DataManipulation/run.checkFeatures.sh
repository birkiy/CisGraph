

strTeth=/groups/lackgrp/ll_members/tunc/phd/ana-starrseq-lncap-lacklab/analysis/bigwig/lncap-logFC-etoh-plasmid.bigWig
strTdht=/groups/lackgrp/ll_members/tunc/phd/ana-starrseq-lncap-lacklab/analysis/bigwig/lncap-logFC-dht-plasmid.bigWig

creDhtData=/home/ualtintas/github/Data/CisGraph/Vers2.0/DHS/DHT_peaks_tabbed.bed
creEthData=/home/ualtintas/github/Data/CisGraph/Vers2.0/DHS/EtOH_peaks_tabbed.bed
feaData=/home/ualtintas/github/Data/CisGraph/Vers2.0/Features
heatmap=/home/ualtintas/github/Data/CisGraph/Vers2.0/Heatmaps

starrseq=$feaData/StarrSeq
groseq=$feaData/GroSeq
cage=$feaData/Cage

strEth=$starrseq/starrseq.etoh.BPM.bigWig
strDht=$starrseq/starrseq.dht.BPM.bigWig


groDmsMin=$groseq/groseq.dmso.BPM.minus.bigWig
groDmsPlu=$groseq/groseq.dmso.BPM.plus.bigWig
groDhtMin=$groseq/groseq.dht.BPM.minus.bigWig
groDhtPlu=$groseq/groseq.dht.BPM.plus.bigWig

cagEthPlu=$cage/cage.+.etoh.BPM.bigWig
cagEthMin=$cage/cage.-.etoh.BPM.bigWig
cagDhtPlu=$cage/cage.+.dht.BPM.bigWig
cagDhtMin=$cage/cage.-.dht.BPM.bigWig

computeMatrix reference-point -S $strDht $strEth -R $creDhtData $creEthData --skipZeros --referencePoint center -a 2000 -b 2000 -o $heatmap/creStarrMat.mat -p 30 -q

plotHeatmap -m $heatmap/creStarrMat.mat -out $heatmap/cerStarrMat.pdf --heatmapWidth 10 --heatmapHeight 50 --colorMap "YlGnBu"




computeMatrix reference-point -S $strDht $strEth $strTeth $strTdht -R $creDhtData $creEthData --skipZeros --referencePoint center -a 2000 -b 2000 -o $heatmap/creStarrMat2.mat -p 30 -q

plotHeatmap -m $heatmap/creStarrMat2.mat -out $heatmap/cerStarrMat2.pdf --heatmapWidth 10 --heatmapHeight 50 --colorMap "YlGnBu"






computeMatrix reference-point -S $groDmsMin $groDmsPlu $groDhtMin $groDhtPlu -R $creDhtData $creEthData --skipZeros --referencePoint center -a 2000 -b 2000 -o $heatmap/creGroMat.mat -p 30 -q

plotHeatmap -m $heatmap/creGroMat.mat -out $heatmap/creGroMat.pdf --heatmapWidth 10 --heatmapHeight 50 --colorMap "YlGnBu"



computeMatrix reference-point -S $cagEthPlu $cagEthMin $cagDhtPlu $cagDhtMin -R $creDhtData $creEthData --skipZeros --referencePoint center -a 2000 -b 2000 -o $heatmap/creCageMat.mat -p 30 -q

plotHeatmap -m $heatmap/creCageMat.mat -out $heatmap/creCageMat.pdf --heatmapWidth 10 --heatmapHeight 50 --colorMap "YlGnBu"




computeMatrix reference-point -S $cagEthPlu $cagEthMin $cagDhtPlu $cagDhtMin $strDht $strEth $groDmsMin $groDmsPlu $groDhtMin $groDhtPlu -R $creDhtData $creEthData --skipZeros --referencePoint center -a 2000 -b 2000 -o $heatmap/creMat.mat -p 30 -q

plotHeatmap -m $heatmap/creMat.mat -out $heatmap/cerMat.pdf --heatmapWidth 10 --heatmapHeight 50 --colorMap "YlGnBu"
