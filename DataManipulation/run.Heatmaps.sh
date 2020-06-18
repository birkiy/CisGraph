#!/bin/bash
#SBATCH --job-name=vers2.bw
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=15
#SBATCH --export=all
#SBATCH -p long


outData=/home/ualtintas/github/Data/CisGraph/Vers2.0/Features/Cage/bigWig
heatmap=/home/ualtintas/github/Data/CisGraph/Vers2.0/Heatmaps
cagePwd=/home/ualtintas/github/Data/CisGraph/Vers2.0/Features/Cage

tsP=/home/ualtintas/genomeAnnotations/Regions/TSS.hg19.+.bed
tsM=/home/ualtintas/genomeAnnotations/Regions/TSS.hg19.-.bed

creDhtData=/home/ualtintas/github/Data/CisGraph/Vers2.0/DHS/DHT_peaks_tabbed.bed
creEthData=/home/ualtintas/github/Data/CisGraph/Vers2.0/DHS/EtOH_peaks_tabbed.bed


groDmsMin=/home/ualtintas/github/Data/CisGraph/Vers2.0/Features/GroSeq/groseq.dmso.-.BPM.bigWig
groDmsPlu=/home/ualtintas/github/Data/CisGraph/Vers2.0/Features/GroSeq/groseq.dmso.+.BPM.bigWig
groDhtMin=/home/ualtintas/github/Data/CisGraph/Vers2.0/Features/GroSeq/groseq.dht.-.BPM.bigWig
groDhtPlu=/home/ualtintas/github/Data/CisGraph/Vers2.0/Features/GroSeq/groseq.dht.+.BPM.bigWig


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

plotProfile --numPlotsPerRow 4 -m $heatmap/tssCageMat.mat -out $heatmap/tssCageMat.pdf --color "#fe8a71" "#3da4ab"

echo "TSS Cage Plot Profile"

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

plotProfile --numPlotsPerRow 4 -m $heatmap/creCageMat5.mat -out $heatmap/creCageMat5.pdf --color "#fb7813" "#17706e"

echo "CRE Cage Plot Profile"

tail -n +4 $cagePwd/creCageMat5.tab > $cagePwd/headles.tab
sed -n '3p' $cagePwd/creCageMat5.tab | cut -f 3- > $cagePwd/header
cat $cagePwd/header $cagePwd/headles.tab > $cagePwd/creCageMat.modified.tab
paste $cagePwd/creCageMat5.bed $cagePwd/creCageMat.modified.tab > $cagePwd/creCageConcat.tab
echo "Concat"

python creCenter.py
# << 'runAfter'

# Run creCenter.py script to align regions on TSS center'

echo "TSS centric CRE"
computeMatrix reference-point -S $groDmsMin $groDmsPlu $groDhtMin $groDhtPlu -R $cagePwd/creCtr.bed $cagePwd/creCtr2.bed $cagePwd/creCtrP.bed $cagePwd/creCtrM.bed --skipZeros --referencePoint center -a 2000 -b 2000 -o $heatmap/creGroCMat.mat -p 30 -q

plotHeatmap -m $heatmap/creGroCMat.mat -out $heatmap/creGroCHMat.pdf --colorMap "YlGnBu" --whatToShow "heatmap and colorbar"

plotProfile --numPlotsPerRow 2 -m $heatmap/creGroCMat.mat -out $heatmap/creGroCPMat.pdf --color "#f6cd61" "#b1b493" "#fe8a71" "#3da4ab"

echo "Gro heatmaps of CREs Plot Profi"

computeMatrix reference-point -S $groDmsMin $groDmsPlu $groDhtMin $groDhtPlu -R $cagePwd/creCtrP.bed $cagePwd/creCtrM.bed --skipZeros --referencePoint center -a 2000 -b 2000 -o $heatmap/creGroCUMat.mat -p 30 -q

plotProfile --numPlotsPerRow 2 -m $heatmap/creGroCUMat.mat -out $heatmap/creGroCUPMat.pdf --color "#fe8a71" "#3da4ab"

# runAfter


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
-R $cagePwd/creCtr.bed $cagePwd/creCtr2.bed $cagePwd/creCtrP.bed $cagePwd/creCtrM.bed --skipZeros --referencePoint center -a 500 -b 500 -o $heatmap/creCageAMat5.mat -p 30 -q

plotProfile --numPlotsPerRow 4 -m $heatmap/creCageAMat5.mat -out $heatmap/creCageAMat5.pdf --color "#f6cd61" "#b1b493" "#fe8a71" "#3da4ab"


hg19=/home/ualtintas/genomeAnnotations/hg19.fa

dht=/home/ualtintas/github/Data/CisGraph/Vers2.0/NodeBeds/creDHT.bed
eth=/home/ualtintas/github/Data/CisGraph/Vers2.0/NodeBeds/creEtOH.bed

dhtS=/home/ualtintas/github/Data/CisGraph/Vers2.0/Sequences/cre.DHT.Seq.fasta
ethS=/home/ualtintas/github/Data/CisGraph/Vers2.0/Sequences/cre.EtOH.Seq.fasta


bedtools getfasta -fi $hg19 -bed $dht -name > $dhtS
bedtools getfasta -fi $hg19 -bed $eth -name > $ethS
