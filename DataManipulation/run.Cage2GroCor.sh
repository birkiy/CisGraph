#!/bin/bash
#SBATCH --job-name=vers2.bw
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=15
#SBATCH --export=all
#SBATCH -p long


outData=/home/ualtintas/github/Data/CisGraph/Vers2.0



multiBigwigSummary bins -b \
$outData/Features/GroSeq/groseq.dht.BPM.minus.bigWig \
$outData/Features/GroSeq/groseq.dht.BPM.plus.bigWig \
$outData/Features/GroSeq/groseq.dmso.BPM.minus.bigWig \
$outData/Features/GroSeq/groseq.dmso.BPM.plus.bigWig \
$outData/Features/Cage/LNCaP.0h.rep1.-.BPM.bw \
$outData/Features/Cage/LNCaP.0h.rep1.+.BPM.bw \
$outData/Features/Cage/LNCaP.3h.rep1.-.BPM.bw \
$outData/Features/Cage/LNCaP.3h.rep1.+.BPM.bw \
$outData/Features/Cage/LNCaP.6h.rep1.-.BPM.bw \
$outData/Features/Cage/LNCaP.6h.rep1.+.BPM.bw \
$outData/Features/Cage/LNCaP.18h.rep1.-.BPM.bw \
$outData/Features/Cage/LNCaP.18h.rep1.+.BPM.bw \
$outData/Features/Cage/LNCaP.48h.rep1.-.BPM.bw \
$outData/Features/Cage/LNCaP.48h.rep1.+.BPM.bw \
$outData/Features/Cage/LNCaP.72h.rep1.-.BPM.bw \
$outData/Features/Cage/LNCaP.72h.rep1.+.BPM.bw \
-o $outData/Cage2GroCor.npz \
--smartLabels --outRawCounts $outData/Cage2GroCor.tab



multiBigwigSummary bins -b \
$outData/Features/GroSeq/groseq.dht.BPM.minus.bigWig \
$outData/Features/GroSeq/groseq.dht.BPM.plus.bigWig \
$outData/Features/GroSeq/groseq.dmso.BPM.minus.bigWig \
$outData/Features/GroSeq/groseq.dmso.BPM.plus.bigWig \
$outData/Features/Cage/LNCaP.0h.rep1.-.BPM.bw \
$outData/Features/Cage/LNCaP.0h.rep1.+.BPM.bw \
$outData/Features/Cage/LNCaP.3h.rep1.-.BPM.bw \
$outData/Features/Cage/LNCaP.3h.rep1.+.BPM.bw \
$outData/Features/Cage/LNCaP.6h.rep1.-.BPM.bw \
$outData/Features/Cage/LNCaP.6h.rep1.+.BPM.bw \
$outData/Features/Cage/LNCaP.18h.rep1.-.BPM.bw \
$outData/Features/Cage/LNCaP.18h.rep1.+.BPM.bw \
$outData/Features/Cage/LNCaP.48h.rep1.-.BPM.bw \
$outData/Features/Cage/LNCaP.48h.rep1.+.BPM.bw \
$outData/Features/Cage/LNCaP.72h.rep1.-.BPM.bw \
$outData/Features/Cage/LNCaP.72h.rep1.+.BPM.bw \
-o $outData/Cage2GroCor.npz \
--smartLabels --outRawCounts $outData/Cage2GroCor.tab




tsP=/home/ualtintas/genomeAnnotations/Regions/TSS.hg19.+.bed
tsM=/home/ualtintas/genomeAnnotations/Regions/TSS.hg19.-.bed

multiBigwigSummary BED-file -b \
$outData/Features/GroSeq/groseq.dht.BPM.minus.bigWig \
$outData/Features/GroSeq/groseq.dht.BPM.plus.bigWig \
$outData/Features/GroSeq/groseq.dmso.BPM.minus.bigWig \
$outData/Features/GroSeq/groseq.dmso.BPM.plus.bigWig \
$outData/Features/Cage/LNCaP.0h.rep1.-.BPM.bw \
$outData/Features/Cage/LNCaP.0h.rep1.+.BPM.bw \
$outData/Features/Cage/LNCaP.3h.rep1.-.BPM.bw \
$outData/Features/Cage/LNCaP.3h.rep1.+.BPM.bw \
$outData/Features/Cage/LNCaP.6h.rep1.-.BPM.bw \
$outData/Features/Cage/LNCaP.6h.rep1.+.BPM.bw \
$outData/Features/Cage/LNCaP.18h.rep1.-.BPM.bw \
$outData/Features/Cage/LNCaP.18h.rep1.+.BPM.bw \
$outData/Features/Cage/LNCaP.48h.rep1.-.BPM.bw \
$outData/Features/Cage/LNCaP.48h.rep1.+.BPM.bw \
$outData/Features/Cage/LNCaP.72h.rep1.-.BPM.bw \
$outData/Features/Cage/LNCaP.72h.rep1.+.BPM.bw \
-o $outData/Cage2GroTsPCor.npz \
--BED $tsP \
--outRawCounts $outData/Cage2GroTsPCor.tab \
-p 30

multiBigwigSummary BED-file -b \
$outData/Features/GroSeq/groseq.dht.BPM.minus.bigWig \
$outData/Features/GroSeq/groseq.dht.BPM.plus.bigWig \
$outData/Features/GroSeq/groseq.dmso.BPM.minus.bigWig \
$outData/Features/GroSeq/groseq.dmso.BPM.plus.bigWig \
$outData/Features/Cage/LNCaP.0h.rep1.-.BPM.bw \
$outData/Features/Cage/LNCaP.0h.rep1.+.BPM.bw \
$outData/Features/Cage/LNCaP.3h.rep1.-.BPM.bw \
$outData/Features/Cage/LNCaP.3h.rep1.+.BPM.bw \
$outData/Features/Cage/LNCaP.6h.rep1.-.BPM.bw \
$outData/Features/Cage/LNCaP.6h.rep1.+.BPM.bw \
$outData/Features/Cage/LNCaP.18h.rep1.-.BPM.bw \
$outData/Features/Cage/LNCaP.18h.rep1.+.BPM.bw \
$outData/Features/Cage/LNCaP.48h.rep1.-.BPM.bw \
$outData/Features/Cage/LNCaP.48h.rep1.+.BPM.bw \
$outData/Features/Cage/LNCaP.72h.rep1.-.BPM.bw \
$outData/Features/Cage/LNCaP.72h.rep1.+.BPM.bw \
-o $outData/Cage2GroTsMCor.npz \
--BED $tsM \
--outRawCounts $outData/Cage2GroTsMCor.tab \
-p 30



plotCorrelation \
   -in Cage2GroCor.npz \
   --corMethod spearman --skipZeros \
   --plotTitle "Spearman Correlation of Read Counts" \
   --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
   -o Cage2GroCor.pdf




plotCorrelation \
  -in Cage2GroTsPCor.npz \
  --corMethod spearman --skipZeros \
  --plotTitle "Spearman Correlation of Read Counts" \
  --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
  -o Cage2GroTsPCor.pdf



plotCorrelation \
 -in Cage2GroTsMCor.npz \
 --corMethod spearman --skipZeros \
 --plotTitle "Spearman Correlation of Read Counts" \
 --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
 -o Cage2GroTsMCor.pdf
