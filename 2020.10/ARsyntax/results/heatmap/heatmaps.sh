#!/bin/bash
#SBATCH --job-name=bpnet.heatmap
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=20
#SBATCH --export=all
#SBATCH -p long


home=/groups/lackgrp/ll_members/berkay/ARsyntax
results=/groups/lackgrp/ll_members/berkay/ARsyntax/results

computeMatrix reference-point -S \
  $home/results/bigwig/22RV1.AR-C19.-.5end.bigWig \
  $home/results/bigwig/22RV1.AR-C19.+.5end.bigWig \
  -R $home/results/peak/idr/22RV1.AR-C19.idr.peaks.bed \
  -a 1000 -b 1000 --sortRegions descend -p 20 --referencePoint=center\
  -o $results/coverage/22RV1.AR-C19.npz

computeMatrix reference-point -S \
  $home/results/bigwig/22RV1.AR-V7.-.5end.bigWig \
  $home/results/bigwig/22RV1.AR-V7.+.5end.bigWig \
  -R $home/results/peak/idr/22RV1.AR-V7.idr.peaks.bed \
  -a 1000 -b 1000 --sortRegions descend -p 20 --referencePoint=center \
  -o $results/coverage/22RV1.AR-V7.npz



plotHeatmap -m $results/coverage/22RV1.AR-C19.npz \
  -out $results/heatmap/22RV1.AR-C19.pdf \
  --whatToShow 'heatmap and colorbar' --colorMap "Blues" --missingDataColor 1 --sortRegions descend \
  --heatmapHeight 20 --heatmapWidth 7

plotHeatmap -m $results/coverage/22RV1.AR-V7.npz \
  -out $results/heatmap/22RV1.AR-V7.pdf \
  --whatToShow 'heatmap and colorbar' --colorMap "Blues" --missingDataColor 1 --sortRegions descend \
  --heatmapHeight 20 --heatmapWidth 7



computeMatrix reference-point -S \
  $home/results/bigwig/LN95.AR-C19.-.5end.bigWig \
  $home/results/bigwig/LN95.AR-C19.+.5end.bigWig \
  -R $home/results/peak/idr/LN95.AR-C19.idr.peaks.bed \
  -a 1000 -b 1000 --sortRegions descend -p 20 --referencePoint=center \
  -o $results/coverage/LN95.AR-C19.npz

computeMatrix reference-point -S \
  $home/results/bigwig/LN95.AR-V7.-.5end.bigWig \
  $home/results/bigwig/LN95.AR-V7.+.5end.bigWig \
  -R $home/results/peak/idr/LN95.AR-V7.idr.peaks.bed \
  -a 1000 -b 1000 --sortRegions descend -p 20 --referencePoint=center \
  -o $results/coverage/LN95.AR-V7.npz

computeMatrix reference-point -S \
  $home/results/bigwig/LNCaP.dht.AR.-.5end.bigWig \
  $home/results/bigwig/LNCaP.dht.AR.+.5end.bigWig \
  -R $home/results/peak/idr/LNCaP.dht.AR.idr.peaks.bed \
  -a 1000 -b 1000 --sortRegions descend -p 20 --referencePoint=center\
  -o $results/coverage/LNCaP.dht.AR.npz

computeMatrix reference-point -S \
  $home/results/bigwig/LNCaP.veh.AR.-.5end.bigWig \
  $home/results/bigwig/LNCaP.veh.AR.+.5end.bigWig \
  -R $home/results/peak/idr/LNCaP.veh.AR.idr.peaks.bed \
  -a 1000 -b 1000 --sortRegions descend -p 20 --referencePoint=center\
  -o $results/coverage/LNCaP.veh.AR.npz

computeMatrix reference-point -S \
  $home/results/bigwig/malignant.1.AR.-.5end.bigWig \
  $home/results/bigwig/malignant.1.AR.+.5end.bigWig \
  -R $home/results/peak/idr/malignant.1.AR.idr.peaks.bed \
  -a 1000 -b 1000 --sortRegions descend -p 20 --referencePoint=center\
  -o $results/coverage/malignant.1.AR.npz

computeMatrix reference-point -S \
  $home/results/bigwig/malignant.2.AR.-.5end.bigWig \
  $home/results/bigwig/malignant.2.AR.+.5end.bigWig \
  -R $home/results/peak/idr/malignant.2.AR.idr.peaks.bed \
  -a 1000 -b 1000 --sortRegions descend -p 20 --referencePoint=center\
  -o $results/coverage/malignant.2.AR.npz

computeMatrix reference-point -S \
  $home/results/bigwig/malignant.3.AR.-.5end.bigWig \
  $home/results/bigwig/malignant.3.AR.+.5end.bigWig \
  -R $home/results/peak/idr/malignant.3.AR.idr.peaks.bed \
  -a 1000 -b 1000 --sortRegions descend -p 20 --referencePoint=center\
  -o $results/coverage/malignant.3.AR.npz

computeMatrix reference-point -S \
  $home/results/bigwig/malignant.4.AR.-.5end.bigWig \
  $home/results/bigwig/malignant.4.AR.+.5end.bigWig \
  -R $home/results/peak/idr/malignant.4.AR.idr.peaks.bed \
  -a 1000 -b 1000 --sortRegions descend -p 20 --referencePoint=center\
  -o $results/coverage/malignant.4.AR.npz

computeMatrix reference-point -S \
  $home/results/bigwig/non-malignant.1.AR.-.5end.bigWig \
  $home/results/bigwig/non-malignant.1.AR.+.5end.bigWig \
  -R $home/results/peak/idr/non-malignant.1.AR.idr.peaks.bed \
  -a 1000 -b 1000 --sortRegions descend -p 20 --referencePoint=center\
  -o $results/coverage/non-malignant.1.AR.npz

computeMatrix reference-point -S \
  $home/results/bigwig/non-malignant.2.AR.-.5end.bigWig \
  $home/results/bigwig/non-malignant.2.AR.+.5end.bigWig \
  -R $home/results/peak/idr/non-malignant.2.AR.idr.peaks.bed \
  -a 1000 -b 1000 --sortRegions descend -p 20 --referencePoint=center\
  -o $results/coverage/non-malignant.2.AR.npz





plotHeatmap -m $results/coverage/LN95.AR-C19.npz \
  -out $results/heatmap/LN95.AR-C19.pdf \
  --whatToShow 'heatmap and colorbar' --colorMap "Blues" --missingDataColor 1 --sortRegions descend \
  --heatmapHeight 20 --heatmapWidth 7

plotHeatmap -m $results/coverage/LN95.AR-V7.npz \
  -out $results/heatmap/LN95.AR-V7.pdf \
  --whatToShow 'heatmap and colorbar' --colorMap "Blues" --missingDataColor 1 --sortRegions descend \
  --heatmapHeight 20 --heatmapWidth 7

plotHeatmap -m $results/coverage/LNCaP.dht.AR.npz \
  -out $results/heatmap/LNCaP.dht.AR.pdf \
  --whatToShow 'heatmap and colorbar' --colorMap "Blues" --missingDataColor 1 --sortRegions descend \
  --heatmapHeight 20 --heatmapWidth 7

plotHeatmap -m $results/coverage/LNCaP.veh.AR.npz \
  -out $results/heatmap/LNCaP.veh.AR.pdf \
  --whatToShow 'heatmap and colorbar' --colorMap "Blues" --missingDataColor 1 --sortRegions descend \
  --heatmapHeight 20 --heatmapWidth 7

plotHeatmap -m $results/coverage/malignant.1.AR.npz \
  -out $results/heatmap/malignant.1.AR.pdf \
  --whatToShow 'heatmap and colorbar' --colorMap "Blues" --missingDataColor 1 --sortRegions descend \
  --heatmapHeight 20 --heatmapWidth 7

plotHeatmap -m $results/coverage/malignant.2.AR.npz \
  -out $results/heatmap/malignant.2.AR.pdf \
  --whatToShow 'heatmap and colorbar' --colorMap "Blues" --missingDataColor 1 --sortRegions descend \
  --heatmapHeight 20 --heatmapWidth 7

plotHeatmap -m $results/coverage/malignant.3.AR.npz \
  -out $results/heatmap/malignant.3.AR.pdf \
  --whatToShow 'heatmap and colorbar' --colorMap "Blues" --missingDataColor 1 --sortRegions descend \
  --heatmapHeight 20 --heatmapWidth 7

plotHeatmap -m $results/coverage/malignant.4.AR.npz \
  -out $results/heatmap/malignant.4.AR.pdf \
  --whatToShow 'heatmap and colorbar' --colorMap "Blues" --missingDataColor 1 --sortRegions descend \
  --heatmapHeight 20 --heatmapWidth 7

plotHeatmap -m $results/coverage/non-malignant.1.AR.npz \
  -out $results/heatmap/non-malignant.1.AR.pdf \
  --whatToShow 'heatmap and colorbar' --colorMap "Blues" --missingDataColor 1 --sortRegions descend \
  --heatmapHeight 20 --heatmapWidth 7

plotHeatmap -m $results/coverage/non-malignant.2.AR.npz \
  -out $results/heatmap/non-malignant.2.AR.pdf \
  --whatToShow 'heatmap and colorbar' --colorMap "Blues" --missingDataColor 1 --sortRegions descend \
  --heatmapHeight 20 --heatmapWidth 7


22RV1.AR-C19.rep1.-.5end.bigWig
22RV1.AR-C19.rep1.+.5end.bigWig
22RV1.AR-C19.rep2.-.5end.bigWig
22RV1.AR-C19.rep2.+.5end.bigWig

22RV1.AR-V7.rep1.-.5end.bigWig
22RV1.AR-V7.rep1.+.5end.bigWig
22RV1.AR-V7.rep2.-.5end.bigWig
22RV1.AR-V7.rep2.+.5end.bigWig

LN95.AR-C19.rep1.-.5end.bigWig
LN95.AR-C19.rep1.+.5end.bigWig
LN95.AR-C19.rep2.-.5end.bigWig
LN95.AR-C19.rep2.+.5end.bigWig

LN95.AR-V7.rep1.-.5end.bigWig
LN95.AR-V7.rep1.+.5end.bigWig
LN95.AR-V7.rep2.-.5end.bigWig
LN95.AR-V7.rep2.+.5end.bigWig

LNCaP.dht.AR.rep1.-.5end.bigWig
LNCaP.dht.AR.rep1.+.5end.bigWig
LNCaP.dht.AR.rep2.-.5end.bigWig
LNCaP.dht.AR.rep2.+.5end.bigWig

LNCaP.veh.AR.rep1.-.5end.bigWig
LNCaP.veh.AR.rep1.+.5end.bigWig
LNCaP.veh.AR.rep2.-.5end.bigWig
LNCaP.veh.AR.rep2.+.5end.bigWig

malignant.1.AR.rep1.-.5end.bigWig
malignant.1.AR.rep1.+.5end.bigWig
malignant.1.AR.rep2.-.5end.bigWig
malignant.1.AR.rep2.+.5end.bigWig

malignant.2.AR.rep1.-.5end.bigWig
malignant.2.AR.rep1.+.5end.bigWig
malignant.2.AR.rep2.-.5end.bigWig
malignant.2.AR.rep2.+.5end.bigWig

malignant.3.AR.rep1.-.5end.bigWig
malignant.3.AR.rep1.+.5end.bigWig
malignant.3.AR.rep2.-.5end.bigWig
malignant.3.AR.rep2.+.5end.bigWig

malignant.4.AR.rep1.-.5end.bigWig
malignant.4.AR.rep1.+.5end.bigWig
malignant.4.AR.rep2.-.5end.bigWig
malignant.4.AR.rep2.+.5end.bigWig

non-malignant.1.AR.rep1.-.5end.bigWig
non-malignant.1.AR.rep1.+.5end.bigWig
non-malignant.1.AR.rep2.-.5end.bigWig
non-malignant.1.AR.rep2.+.5end.bigWig

non-malignant.2.AR.rep1.-.5end.bigWig
non-malignant.2.AR.rep1.+.5end.bigWig
non-malignant.2.AR.rep2.-.5end.bigWig
non-malignant.2.AR.rep2.+.5end.bigWig
