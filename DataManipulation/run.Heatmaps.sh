

outData=/home/ualtintas/github/Data/CisGraph/Vers2.0/Features/Cage/bigWig
heatmap=/home/ualtintas/github/Data/CisGraph/Vers2.0/Heatmaps
cagePwd=/home/ualtintas/github/Data/CisGraph/Vers2.0/Features/Cage

tsP=/home/ualtintas/genomeAnnotations/Regions/TSS.hg19.+.bed
tsM=/home/ualtintas/genomeAnnotations/Regions/TSS.hg19.-.bed

creDhtData=/home/ualtintas/github/Data/CisGraph/Vers2.0/DHS/DHT_peaks_tabbed.bed
creEthData=/home/ualtintas/github/Data/CisGraph/Vers2.0/DHS/EtOH_peaks_tabbed.bed


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

tail -n +4 $cagePwd/creCageMat5.tab > $cagePwd/headles.tab
sed -n '3p' $cagePwd/creCageMat5.tab | cut -f 3- > $cagePwd/header
cat $cagePwd/header $cagePwd/headles.tab > $cagePwd/creCageMat.modified.tab
paste $cagePwd/creCageMat5.bed $cagePwd/creCageMat.modified.tab > $cagePwd/creCageConcat.tab


python creCenter.py
# << 'runAfter'

# Run creCenter.py script to align regions on TSS center'


computeMatrix reference-point -S $groDmsMin $groDmsPlu $groDhtMin $groDhtPlu -R $cage/creCtr.bed $cage/creCtrP.bed $cage/creCtrM.bed --skipZeros --referencePoint center -a 2000 -b 2000 -o $heatmap/creGroCMat.mat -p 30 -q

plotHeatmap -m $heatmap/creGroCMat.mat -out $heatmap/creGroCMat.pdf --heatmapWidth 10 --heatmapHeight 50 --colorMap "YlGnBu"



computeMatrix reference-point -S $groDmsMin $groDmsPlu $groDhtMin $groDhtPlu -R $cage/creCtr.bed --skipZeros --referencePoint center -a 2000 -b 2000 -o $heatmap/creGroCBMat.mat -p 30 -q

plotHeatmap -m $heatmap/creGroCBMat.mat -out $heatmap/creGroCBMat.pdf --heatmapWidth 10 --heatmapHeight 50 --colorMap "YlGnBu"


computeMatrix reference-point -S $groDmsMin $groDmsPlu $groDhtMin $groDhtPlu -R $cage/creCtrP.bed $cage/creCtrM.bed --skipZeros --referencePoint center -a 2000 -b 2000 -o $heatmap/creGroCDMat.mat -p 30 -q

plotHeatmap -m $heatmap/creGroCDMat.mat -out $heatmap/creGroCDMat.pdf --heatmapWidth 10 --heatmapHeight 50 --colorMap "YlGnBu"

# runAfter
