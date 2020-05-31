


mi0=~/BigWig/miRNA.0h.DHT.bigWig
mi6=~/BigWig/miRNA.6h.DHT.bigWig
mi24=~/BigWig/miRNA.24h.DHT.bigWig

mi0=~/LncapExo/miRNA/miRNA.0h.DHT.bigWig
mi6=~/LncapExo/miRNA/miRNA.6h.DHT.bigWig
mi24=~/LncapExo/miRNA/miRNA.24h.DHT.bigWig


tsP=~/genomeAnnotations/Regions/TSS.hg19.+.bed
tsM=~/genomeAnnotations/Regions/TSS.hg19.-.bed


computeMatrix reference-point -S $mi0 $mi6 $mi24 -R $tsP $tsM --referencePoint center -a 10000 -b 10000 -o mat/miR.TSS.10.mat.gz -p 30

plotHeatmap -m mat/miR.TSS.10.mat.gz -out miR.TSS.10.pdf --heatmapWidth 10 --heatmapHeight 70 --colorMap "YlGnBu"


computeMatrix reference-point -S $mi0 $mi6 $mi24 -R $tsP $tsM --binSize 5000 --referencePoint center -a 500000 -b 500000 -o mat/miR.TSS.500.mat.gz -p 30

plotHeatmap -m mat/miR.TSS.500.mat.gz -out miR.TSS.500.pdf --heatmapWidth 10 --heatmapHeight 70 --colorMap "YlGnBu"




ttP=~/genomeAnnotations/Regions/TTS.hg19.+.bed
ttM=~/genomeAnnotations/Regions/TTS.hg19.-.bed


computeMatrix reference-point -S $mi0 $mi6 $mi24 -R $ttP $ttM --referencePoint center -a 10000 -b 10000 -o mat/miR.TTS.10.mat.gz -p 20

plotHeatmap -m mat/miR.TTS.10.mat.gz -out miR.TTS.10.pdf --heatmapWidth 10 --heatmapHeight 70 --colorMap "YlGnBu"

computeMatrix reference-point -S $mi0 $mi6 $mi24 -R $ttP $ttM --binSize 5000 --referencePoint center -a 500000 -b 500000 -o mat/miR.TTS.500.mat.gz -p 30

plotHeatmap -m mat/miR.TTS.500.mat.gz -out miR.TTS.500.pdf --heatmapWidth 10 --heatmapHeight 70 --colorMap "YlGnBu"



con=~/ARBSs/regions/cons-arbs.bed
ind=~/ARBSs/regions/ind-arbs.bed
non=~/ARBSs/regions/Non-Active-ARBS.bed


computeMatrix reference-point -S $mi0 $mi6 $mi24 -R $con $ind $non --referencePoint center -a 100000 -b 100000 -o mat/miR.ARBS.100.mat.gz -p 30

plotHeatmap -m mat/miR.ARBS.100.mat.gz -out miR.ARBS.100.pdf --heatmapWidth 10 --heatmapHeight 70 --colorMap "YlGnBu"



computeMatrix reference-point -S $mi0 $mi6 $mi24 -R  $con $ind $non --binSize 5000 --referencePoint center -a 1000000 -b 1000000 -o mat/miR.ARBS.1000.mat.gz -p 30

plotHeatmap -m mat/miR.ARBS.1000.mat.gz -out miR.ARBS.1000.pdf --heatmapWidth 10 --heatmapHeight 70 --colorMap "YlGnBu"




computeMatrix reference-point -S $mi0 $mi6 $mi24 -R $con $ind $non --referencePoint center -a 60000 -b 60000 -o mat/miR.ARBS.60.mat.gz -p 20

plotHeatmap -m mat/miR.ARBS.60.mat.gz -out miR.ARBS.60.pdf --heatmapWidth 10 --heatmapHeight 70 --colorMap "YlGnBu"


parts=~/ARBSs/Parts/*

computeMatrix scale-regions -S $mi0 $mi6 $mi24 -R $parts -a 0 -b 0 --binSize 60000 --transcript_id_designator ACTB --outFileNameMatrix PlusMinusCREMat.tab --outFileSortedRegions PlusMinusCREMatRegions.bed -o mat/groCREPartDensity.mat.gz -p 30 -m 350

plotHeatmap -m mat/groCREPartDensity.mat.gz -out groHeatmapCREParts.pdf --heatmapWidth 10 --heatmapHeight 70 --colorMap "YlGnBu"

cd tmp

tail -n +4 ../PlusMinusCREMat.tab > CREMat.tab
echo -e "groseq_dht.plus\tgroseq_dht.minus\tgroseq_dmso.plus\tgroseq_dmso.minus"  | cat - CREMat.tab > headCREmat.tab
paste ../PlusMinusCREMatRegions.bed headCREmat.tab > concatCRE.tab


computeMatrix scale-regions -S $mi0 $mi6 $mi24 -R sorted.ttM.tsM.bed -a 3000 -b 3000 --binSize 50 --transcript_id_designator ACTB --outFileNameMatrix scaleTt.ts.tab --outFileSortedRegions scaleTt.ts.Regions.bed -o mat/scaleTt.ts.mat.gz -p 30 -m 5000 --startLabel TTS --endLabel TSS

plotHeatmap -m mat/scaleTt.ts.mat.gz -out scaleTt.ts.pdf --heatmapWidth 10 --heatmapHeight 70 --colorMap "YlGnBu"
