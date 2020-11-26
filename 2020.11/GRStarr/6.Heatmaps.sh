#!/bin/bash
#SBATCH --job-name=starr.heatmaps
#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task 20
#SBATCH --export=all
#SBATCH -p express


LFC=(
"0.5"
"0.6"
"0.7"
"0.8"
"0.9"
"1"
)

DUR=(
"1h"
"4h"
"8h"
"12h"
)

results=/groups/lackgrp/ll_members/berkay/STARRbegin

for lfc in {0..5}
do
  for dur in {0..3}
  do
    # con=$results/peaks/con.${LFC[$lfc]}.GR.${DUR[$dur]}.bed;
    # ind=$results/peaks/ind.${LFC[$lfc]}.GR.${DUR[$dur]}.bed;
    # non=$results/peaks/non.${LFC[$lfc]}.GR.${DUR[$dur]}.bed;
    # nAR=$results/peaks/negativeControlGR.final.bed;
    # computeMatrix reference-point -S \
    #   $results/results/bigwig/GR.0h.dex.merged.RPKM.bigWig \
    #   $results/results/bigwig/GR.1h.dex.merged.RPKM.bigWig \
    #   $results/results/bigwig/GR.4h.dex.merged.RPKM.bigWig \
    #   $results/results/bigwig/GR.8h.dex.merged.RPKM.bigWig \
    #   $results/results/bigwig/GR.12h.dex.merged.RPKM.bigWig \
    #   -R \
    #     $con \
    #     $ind \
    #     $non \
    #     $nAR \
    #   -a 1000 -b 1000 --sortRegions descend -p 20 --referencePoint=center\
    #   -o $results/results/coverage/mat.${LFC[$lfc]}.${DUR[$dur]}.npz;
    plotHeatmap -m $results/results/coverage/mat.${LFC[$lfc]}.${DUR[$dur]}.npz \
      --colorMap "Blues" --missingDataColor 1 \
      --sortRegions descend \
      -out $results/results/coverage/mat.${LFC[$lfc]}.${DUR[$dur]}.pdf \
      --samplesLabel "0h" "1h" "4h" "8h" "12h" --xAxisLabel "Distance(kb)" \
      --heatmapHeight 20 --heatmapWidth 7
  done
done
