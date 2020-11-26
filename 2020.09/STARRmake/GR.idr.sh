



home=/groups/lackgrp/ll_members/berkay/STARRbegin/results/peaks

gunzip $home/GSM2095218_DEX_3hr_GR_ChIPseq.Rep1_peaks.bed.gz
gunzip $home/GSM2095219_DEX_3hr_GR_ChIPseq.Rep2_peaks.bed.gz



sort -k5,5nr $home/GSM2095218_DEX_3hr_GR_ChIPseq.Rep1_peaks.bed > $home/GSM2095218_DEX_3hr_GR_ChIPseq.Rep1_peaks.sorted.bed
sort -k5,5nr $home/GSM2095219_DEX_3hr_GR_ChIPseq.Rep2_peaks.bed > $home/GSM2095219_DEX_3hr_GR_ChIPseq.Rep2_peaks.sorted.bed

idr --samples $home/GSM2095218_DEX_3hr_GR_ChIPseq.Rep1_peaks.sorted.bed $home/GSM2095219_DEX_3hr_GR_ChIPseq.Rep2_peaks.sorted.bed \
--input-file-type bed \
--rank 5 \
--output-file $home/GR.dex.idr \
--plot



bedtools intersect -a $home/GSM2095218_DEX_3hr_GR_ChIPseq.Rep1_peaks.bed -b $home/GSM2095219_DEX_3hr_GR_ChIPseq.Rep2_peaks.bed -u > $home/common.GR.peaks.bed
cut -f1,2,3,4 $home/common.GR.peaks.bed > $home/common.cut.GR.peaks.bed

bedtools intersect -a $home/GSM2095218_DEX_3hr_GR_ChIPseq.Rep1_peaks.bed -b $home/GSM2095219_DEX_3hr_GR_ChIPseq.Rep2_peaks.bed | grep "MACS_peak_124"








bamCoverage --bam \
  $home/mapping/processed/GR.0h.merged.bam \
  -o $home/bigwig/GR.0h.merged.bigWig \
  --extendReads 150 -p 20\
  --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed

bamCoverage --bam \
  $home/mapping/processed/GR.1h.merged.bam \
  -o $home/bigwig/GR.1h.merged.bigWig \
  --extendReads 150 -p 20\
  --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed

bamCoverage --bam \
  $home/mapping/processed/GR.4h.merged.bam \
  -o $home/bigwig/GR.4h.merged.bigWig \
  --extendReads 150 -p 20\
  --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed

bamCoverage --bam \
  $home/mapping/processed/GR.8h.merged.bam \
  -o $home/bigwig/GR.8h.merged.bigWig \
  --extendReads 150 -p 20\
  --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed

bamCoverage --bam \
  $home/mapping/processed/GR.12h.merged.bam \
  -o $home/bigwig/GR.12h.merged.bigWig \
  --extendReads 150 -p 20\
  --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed

bamCoverage --bam \
$home/mapping/processed/input.GSE114063.merged.bam \
-o $home/bigwig/input.GSE114063.merged.bigWig \
--extendReads 150 -p 20\
--blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed





bigWigOverbed
computeMatrix reference-point -S \
  $home/bigvig/GR.0h.merged.bigWig \
  $home/bigvig/GR.4h.merged.bigWig \
  $home/bigvig/GR.8h.merged.bigWig \
  $home/bigvig/GR.8h.merged.bigWig \
  $home/mbigvigGR.12h.merged.bigWig \
  -R <bed-file> \
  -a 1000 -b 1000 --sortRegions descend -p 20 --referencePoint=center\
  -o $results/coverage/matRPKM.npz

$home/mapping/processed/GR.12h.dex.rep1.final.bam


bamCoverage --bam $results/mapping/processed/GR.0h.dex.rep1.final.bam  -o $results/bigwig/GR.0h.dex.rep1.extended.RPKM.bigWig \
  --extendReads 150 \
  --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed \
  -p 20


bamCoverage --bam $results/mapping/processed/GR.12h.dex.rep1.final.bam \
  -o $results/bigwig/GR.12h.dex.rep1.extended.RPKM.bigWig \
  --extendReads 150 \
  --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed \
  -p 20

computeMatrix reference-point -S \
  $results/bigwig/GR.0h.dex.rep1.extended.RPKM.bigWig \
  $results/bigwig/GR.12h.dex.rep1.extended.RPKM.bigWig \
  -R \
  $results/peaks/common.GR.peaks.bed \
  -a 1000 -b 1000 --sortRegions descend -p 20 --referencePoint=center\
  -o $results/coverage/mat0hRPKM.npz


plotHeatmap -m $results/coverage/mat0hRPKM.npz \
  --whatToShow 'heatmap and colorbar' \
  -out $results/coverage/mat0hRPKM.pdf

results=/groups/lackgrp/ll_members/berkay/STARRbegin/results


computeMatrix reference-point -S \
  $results/bigwig/GR.0h.dex.rep1.extended.forward.RPKM.bigWig \
  $results/bigwig/GR.0h.dex.rep1.extended.reverse.RPKM.bigWig \
  $results/bigwig/GR.0h.dex.rep2.extended.forward.RPKM.bigWig \
  $results/bigwig/GR.0h.dex.rep2.extended.reverse.RPKM.bigWig \
  $results/bigwig/GR.0h.dex.rep3.extended.forward.RPKM.bigWig \
  $results/bigwig/GR.0h.dex.rep3.extended.reverse.RPKM.bigWig \
  $results/bigwig/GR.0h.dex.rep4.extended.forward.RPKM.bigWig \
  $results/bigwig/GR.0h.dex.rep4.extended.reverse.RPKM.bigWig \
  $results/bigwig/GR.0h.dex.rep5.extended.forward.RPKM.bigWig \
  $results/bigwig/GR.0h.dex.rep5.extended.reverse.RPKM.bigWig \
  -R \
  $results/peaks/common.GR.peaks.bed \
  -a 1000 -b 1000 --sortRegions descend -p 20 \
  -o $results/coverage/mat0hRPKM.npz



# 1h
computeMatrix reference-point -S \
  $results/bigwig/GR.1h.dex.rep1.extended.forward.RPKM.bigWig \
  $results/bigwig/GR.1h.dex.rep1.extended.reverse.RPKM.bigWig \
  $results/bigwig/GR.1h.dex.rep2.extended.forward.RPKM.bigWig \
  $results/bigwig/GR.1h.dex.rep2.extended.reverse.RPKM.bigWig \
  $results/bigwig/GR.1h.dex.rep3.extended.forward.RPKM.bigWig \
  $results/bigwig/GR.1h.dex.rep3.extended.reverse.RPKM.bigWig \
  $results/bigwig/GR.1h.dex.rep4.extended.forward.RPKM.bigWig \
  $results/bigwig/GR.1h.dex.rep4.extended.reverse.RPKM.bigWig \
  $results/bigwig/GR.1h.dex.rep5.extended.forward.RPKM.bigWig \
  $results/bigwig/GR.1h.dex.rep5.extended.reverse.RPKM.bigWig \
  -R \
  $results/coverage/GR/GRdeInd.bed \
  $results/coverage/GR/GRdeCon.bed \
  $results/coverage/GR/GRdeNon.bed \
  -a 1000 -b 1000 --sortRegions descend -p 20 \
  -o $results/coverage/mat1hRPKM.npz

# 4h
computeMatrix reference-point -S \
  $results/bigwig/GR.4h.dex.rep1.extended.forward.RPKM.bigWig \
  $results/bigwig/GR.4h.dex.rep1.extended.reverse.RPKM.bigWig \
  $results/bigwig/GR.4h.dex.rep2.extended.forward.RPKM.bigWig \
  $results/bigwig/GR.4h.dex.rep2.extended.reverse.RPKM.bigWig \
  $results/bigwig/GR.4h.dex.rep3.extended.forward.RPKM.bigWig \
  $results/bigwig/GR.4h.dex.rep3.extended.reverse.RPKM.bigWig \
  $results/bigwig/GR.4h.dex.rep4.extended.forward.RPKM.bigWig \
  $results/bigwig/GR.4h.dex.rep4.extended.reverse.RPKM.bigWig \
  $results/bigwig/GR.4h.dex.rep5.extended.forward.RPKM.bigWig \
  $results/bigwig/GR.4h.dex.rep5.extended.reverse.RPKM.bigWig \
  -R \
  $results/coverage/GR/GRdeInd.bed \
  $results/coverage/GR/GRdeCon.bed \
  $results/coverage/GR/GRdeNon.bed \
  -a 1000 -b 1000 --sortRegions descend -p 20 \
  -o $results/coverage/mat4hRPKM.npz



# 8h
computeMatrix reference-point -S \
  $results/bigwig/GR.8h.dex.rep1.extended.forward.RPKM.bigWig \
  $results/bigwig/GR.8h.dex.rep1.extended.reverse.RPKM.bigWig \
  $results/bigwig/GR.8h.dex.rep2.extended.forward.RPKM.bigWig \
  $results/bigwig/GR.8h.dex.rep2.extended.reverse.RPKM.bigWig \
  $results/bigwig/GR.8h.dex.rep3.extended.forward.RPKM.bigWig \
  $results/bigwig/GR.8h.dex.rep3.extended.reverse.RPKM.bigWig \
  $results/bigwig/GR.8h.dex.rep4.extended.forward.RPKM.bigWig \
  $results/bigwig/GR.8h.dex.rep4.extended.reverse.RPKM.bigWig \
  $results/bigwig/GR.8h.dex.rep5.extended.forward.RPKM.bigWig \
  $results/bigwig/GR.8h.dex.rep5.extended.reverse.RPKM.bigWig \
  -R \
  $results/coverage/GR/GRdeInd.bed \
  $results/coverage/GR/GRdeCon.bed \
  $results/coverage/GR/GRdeNon.bed \
  -a 1000 -b 1000 --sortRegions descend -p 20 \
  -o $results/coverage/mat8hRPKM.npz



# 12h
computeMatrix reference-point -S \
  $results/bigwig/GR.12h.dex.rep1.extended.forward.RPKM.bigWig \
  $results/bigwig/GR.12h.dex.rep1.extended.reverse.RPKM.bigWig \
  $results/bigwig/GR.12h.dex.rep2.extended.forward.RPKM.bigWig \
  $results/bigwig/GR.12h.dex.rep2.extended.reverse.RPKM.bigWig \
  $results/bigwig/GR.12h.dex.rep3.extended.forward.RPKM.bigWig \
  $results/bigwig/GR.12h.dex.rep3.extended.reverse.RPKM.bigWig \
  $results/bigwig/GR.12h.dex.rep4.extended.forward.RPKM.bigWig \
  $results/bigwig/GR.12h.dex.rep4.extended.reverse.RPKM.bigWig \
  $results/bigwig/GR.12h.dex.rep5.extended.forward.RPKM.bigWig \
  $results/bigwig/GR.12h.dex.rep5.extended.reverse.RPKM.bigWig \
  -R \
  $results/peaks/common.GR.peaks.bed \
  -a 1000 -b 1000 --sortRegions descend -p 20 \
  -o $results/coverage/mat12hRPKM.npz


plotHeatmap -m $results/coverage/mat0hRPKM.npz \
  --whatToShow 'heatmap and colorbar' \
  -out $results/coverage/mat0hRPKM.pdf

plotHeatmap -m $results/coverage/mat1hRPKM.npz \
  --whatToShow 'heatmap and colorbar' \
  -out $results/coverage/mat1hRPKM.pdf


plotHeatmap -m $results/coverage/mat4hRPKM.npz \
  --whatToShow 'heatmap and colorbar' \
  -out $results/coverage/mat4hRPKM.pdf

plotHeatmap -m $results/coverage/mat8hRPKM.npz \
  --whatToShow 'heatmap and colorbar' \
  -out $results/coverage/mat8hRPKM.pdf

plotHeatmap -m $results/coverage/mat12hRPKM.npz \
  --whatToShow 'heatmap and colorbar' \
  -out $results/coverage/mat12hRPKM.pdf


# GR.0h.dex.rep1.extended.RPKM.bigWig GR.0h.dex.rep1.extended.RPKM.bigWig
# GR.0h.dex.rep2.extended.RPKM.bigWig GR.0h.dex.rep2.extended.RPKM.bigWig
# GR.0h.dex.rep3.extended.RPKM.bigWig GR.0h.dex.rep3.extended.RPKM.bigWig
# GR.0h.dex.rep4.extended.RPKM.bigWig GR.0h.dex.rep4.extended.RPKM.bigWig
# GR.0h.dex.rep5.extended.RPKM.bigWig GR.0h.dex.rep5.extended.RPKM.bigWig
# GR.1h.dex.rep1.extended.RPKM.bigWig GR.1h.dex.rep1.extended.RPKM.bigWig
# GR.1h.dex.rep2.extended.RPKM.bigWig GR.1h.dex.rep2.extended.RPKM.bigWig
# GR.1h.dex.rep3.extended.RPKM.bigWig GR.1h.dex.rep3.extended.RPKM.bigWig
# GR.1h.dex.rep4.extended.RPKM.bigWig GR.1h.dex.rep4.extended.RPKM.bigWig
# GR.1h.dex.rep5.extended.RPKM.bigWig GR.1h.dex.rep5.extended.RPKM.bigWig
# GR.4h.dex.rep1.extended.RPKM.bigWig GR.4h.dex.rep1.extended.RPKM.bigWig
# GR.4h.dex.rep2.extended.RPKM.bigWig GR.4h.dex.rep2.extended.RPKM.bigWig
# GR.4h.dex.rep3.extended.RPKM.bigWig GR.4h.dex.rep3.extended.RPKM.bigWig
# GR.4h.dex.rep4.extended.RPKM.bigWig GR.4h.dex.rep4.extended.RPKM.bigWig
# GR.4h.dex.rep5.extended.RPKM.bigWig GR.4h.dex.rep5.extended.RPKM.bigWig
# GR.8h.dex.rep1.extended.RPKM.bigWig GR.8h.dex.rep1.extended.RPKM.bigWig
# GR.8h.dex.rep2.extended.RPKM.bigWig GR.8h.dex.rep2.extended.RPKM.bigWig
# GR.8h.dex.rep3.extended.RPKM.bigWig GR.8h.dex.rep3.extended.RPKM.bigWig
# GR.8h.dex.rep4.extended.RPKM.bigWig GR.8h.dex.rep4.extended.RPKM.bigWig
# GR.8h.dex.rep5.extended.RPKM.bigWig GR.8h.dex.rep5.extended.RPKM.bigWig
# GR.12h.dex.rep1.extended.RPKM.bigWig GR.12h.dex.rep1.extended.RPKM.bigWig
# GR.12h.dex.rep2.extended.RPKM.bigWig GR.12h.dex.rep2.extended.RPKM.bigWig
# GR.12h.dex.rep3.extended.RPKM.bigWig GR.12h.dex.rep3.extended.RPKM.bigWig
# GR.12h.dex.rep4.extended.RPKM.bigWig GR.12h.dex.rep4.extended.RPKM.bigWig
# GR.12h.dex.rep5.extended.RPKM.bigWig GR.12h.dex.rep5.extended.RPKM.bigWig



/groups/lackgrp/ll_members/berkay/STARRbegin/results/coverage/mat12hRPKM.pdf
