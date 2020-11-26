


results=/groups/lackgrp/ll_members/berkay/STARRbegin/results

computeMatrix reference-point -S \
  $results/bigwig/GR.0h.dex.plasmid.merged.bigWig \
  $results/bigwig/GR.1h.dex.plasmid.merged.bigWig \
  $results/bigwig/GR.4h.dex.plasmid.merged.bigWig \
  $results/bigwig/GR.8h.dex.plasmid.merged.bigWig \
  $results/bigwig/GR.12h.dex.plasmid.merged.bigWig \
  -R \
  $results/peaks/common.GR.peaks.bed \
  $results/peaks/negativeControlGR.final.bed \
  -a 1000 -b 1000 --sortRegions descend -p 20 --referencePoint=center\
  -o $results/coverage/matPlasmid.npz


plotHeatmap -m $results/coverage/matPlasmid.npz \
  --whatToShow 'heatmap and colorbar' --colorMap "Blues" --missingDataColor 1 \
  --sortRegions descend \
  -out $results/coverage/matPlasmid.pdf \
  --samplesLabel "0h" "1h" "4h" "8h" "12h" --regionsLabel "GR Binding Sites" "Non GR Binding Sites with GRE" --xAxisLabel "Distance(kb)" \
  --heatmapHeight 20 --heatmapWidth 7




#!/bin/bash
#SBATCH --job-name=plot
#SBATCH --time=23:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=20
#SBATCH --export=all
#SBATCH -p long

results=/groups/lackgrp/ll_members/berkay/STARRbegin/results


computeMatrix reference-point -S \
  $results/bigwig/GR.1h.dex.merged.RPKM.SES.bigWig \
  $results/bigwig/GR.4h.dex.merged.RPKM.SES.bigWig \
  $results/bigwig/GR.8h.dex.merged.RPKM.SES.bigWig \
  $results/bigwig/GR.12h.dex.merged.RPKM.SES.bigWig \
  $results/bigwig/GR.0h.dex.merged.RPKM.SES.bigWig \
  -R \
    /groups/lackgrp/ll_members/berkay/STARRbegin/peaks/common.GR.peaks.bed \
    /groups/lackgrp/ll_members/berkay/STARRbegin/peaks/negativeControlGR.final.bed \
  --referencePoint=center\
  -a 1000 -b 1000 \
  --sortRegions descend -p 16 \
  -o $results/coverage/coverage.RPKM.SES.npz

plotHeatmap -m $results/coverage/coverage.RPKM.SES.npz \
--whatToShow 'heatmap and colorbar' \
--colorMap "Blues" --missingDataColor 1 \
--sortRegions descend \
--heatmapHeight 20 --heatmapWidth 7 \
-out $results/coverage/plot.Heatmap.RPKM.SES.pdf

plotProfile -m $results/coverage/coverage.RPKM.SES.npz \
-out $results/coverage/plot.Profile.RPKM.SES.pdf \
--numPlotsPerRow 2


#!/bin/bash
#SBATCH --job-name=plot
#SBATCH --time=23:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=20
#SBATCH --export=all
#SBATCH -p long

results=/groups/lackgrp/ll_members/berkay/STARRbegin/results

computeMatrix reference-point -S \
  $results/bigwig/GR.8h.dex.merged.RPKM.bigWig \
  $results/bigwig/GR.0h.dex.merged.RPKM.bigWig \
  -R \
    /groups/lackgrp/ll_members/berkay/STARRbegin/peaks/con.GR.bed \
    /groups/lackgrp/ll_members/berkay/STARRbegin/peaks/ind.GR.bed \
    /groups/lackgrp/ll_members/berkay/STARRbegin/peaks/non.GR.bed \
  --referencePoint=center\
  -a 1000 -b 1000 \
  --sortRegions descend -p 16 \
  -o $results/coverage/coverage.nodeClass.npz

plotHeatmap -m $results/coverage/coverage.nodeClass.npz\
  --whatToShow 'heatmap and colorbar' \
  --colorMap "Blues" --missingDataColor 1 \
  --sortRegions descend \
  --heatmapHeight 20 --heatmapWidth 7 \
  --xAxisLabel "Distance(kb)" \
  -out $results/coverage/plot.Heatmap.nodeClass.pdf

plotProfile -m $results/coverage/coverage.nodeClass.npz \
-out $results/coverage/plot.Profile.nodeClass.pdf \
--numPlotsPerRow 2
