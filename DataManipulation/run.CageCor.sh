

outData=/home/ualtintas/github/Data/CisGraph/Vers2.0/Features/Cage


tsP=/home/ualtintas/genomeAnnotations/Regions/TSS.hg19.+.bed
tsM=/home/ualtintas/genomeAnnotations/Regions/TSS.hg19.-.bed

multiBigwigSummary BED-file -b \
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
-o $outData/Cage2GroTsPCor.npz \
--BED $tsP \
--outRawCounts $outData/Cage2GroTsPCor.tab \
-p 30



plotCorrelation \
-in Cage2GroTsPCor.npz \
--corMethod spearman --skipZeros \
--plotTitle "Spearman Correlation of Read Counts" \
--whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
-o Cage2GroTsPCor.pdf



multiBigwigSummary BED-file -b \
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
-o $outData/Cage2GroTsMCor.npz \
--BED $tsM \
--outRawCounts $outData/Cage2GroTsMCor.tab \
-p 30


plotCorrelation \
-in Cage2GroTsMCor.npz \
--corMethod spearman --skipZeros \
--plotTitle "Spearman Correlation of Read Counts" \
--whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
-o Cage2GroTsMCor.pdf
