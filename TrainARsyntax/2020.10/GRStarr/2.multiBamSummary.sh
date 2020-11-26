#!/bin/bash
#SBATCH --job-name=starr.count
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task 64
#SBATCH --export=all
#SBATCH -p long


# home=/groups/lackgrp/ll_members/berkay/STARRbegin/results/
#
# multiBamSummary BED-file \
#   --BED $home/peaks/totalGRE.bed --bamfiles \
#   $home/mapping/processed/GR.0h.dex.rep1.final.bam \
#   $home/mapping/processed/GR.0h.dex.rep2.final.bam \
#   $home/mapping/processed/GR.0h.dex.rep3.final.bam \
#   $home/mapping/processed/GR.0h.dex.rep4.final.bam \
#   $home/mapping/processed/GR.0h.dex.rep5.final.bam \
#   $home/mapping/processed/GR.12h.dex.rep1.final.bam \
#   $home/mapping/processed/GR.12h.dex.rep2.final.bam \
#   $home/mapping/processed/GR.12h.dex.rep3.final.bam \
#   $home/mapping/processed/GR.12h.dex.rep4.final.bam \
#   $home/mapping/processed/GR.12h.dex.rep5.final.bam \
#   $home/mapping/processed/GR.1h.dex.rep1.final.bam \
#   $home/mapping/processed/GR.1h.dex.rep2.final.bam \
#   $home/mapping/processed/GR.1h.dex.rep3.final.bam \
#   $home/mapping/processed/GR.1h.dex.rep4.final.bam \
#   $home/mapping/processed/GR.1h.dex.rep5.final.bam \
#   $home/mapping/processed/GR.4h.dex.rep1.final.bam \
#   $home/mapping/processed/GR.4h.dex.rep2.final.bam \
#   $home/mapping/processed/GR.4h.dex.rep3.final.bam \
#   $home/mapping/processed/GR.4h.dex.rep4.final.bam \
#   $home/mapping/processed/GR.4h.dex.rep5.final.bam \
#   $home/mapping/processed/GR.8h.dex.rep1.final.bam \
#   $home/mapping/processed/GR.8h.dex.rep2.final.bam \
#   $home/mapping/processed/GR.8h.dex.rep3.final.bam \
#   $home/mapping/processed/GR.8h.dex.rep4.final.bam \
#   $home/mapping/processed/GR.8h.dex.rep5.final.bam \
#   $home/mapping/processed/input.GSE114063.pool10.final.bam \
#   $home/mapping/processed/input.GSE114063.pool11.final.bam \
#   $home/mapping/processed/input.GSE114063.pool12.final.bam \
#   $home/mapping/processed/input.GSE114063.pool1.final.bam \
#   $home/mapping/processed/input.GSE114063.pool2.final.bam \
#   $home/mapping/processed/input.GSE114063.pool3.final.bam \
#   $home/mapping/processed/input.GSE114063.pool4.final.bam \
#   $home/mapping/processed/input.GSE114063.pool5.final.bam \
#   $home/mapping/processed/input.GSE114063.pool6.final.bam \
#   $home/mapping/processed/input.GSE114063.pool7.final.bam \
#   $home/mapping/processed/input.GSE114063.pool8.final.bam \
#   $home/mapping/processed/input.GSE114063.pool9.final.bam \
#   --extendReads --samFlagInclude 64 \
#   --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed \
#   -p 64 --outRawCounts totalCountTable2.tsv -out totalCountTable2.npz




home=/groups/lackgrp/ll_members/berkay/STARRbegin

multiBamSummary BED-file \
  --BED /groups/lackgrp/ll_members/berkay/STARRbegin/peaks/totalGRE.bed \
  --bamfiles \
    $home/results/mapping/processed/GR.0h.dex.merged.final.bam \
    $home/results/mapping/processed/GR.12h.dex.merged.final.bam \
    $home/results/mapping/processed/GR.1h.dex.merged.final.bam \
    $home/results/mapping/processed/GR.4h.dex.merged.final.bam \
    $home/results/mapping/processed/GR.8h.dex.merged.final.bam \
    $home/results/mapping/processed/input.GSE114063.merged.final.bam \
  --extendReads --samFlagInclude 64 \
  --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed \
  --centerReads \
  -p 64 \
  --outRawCounts $home/results/coverage/countMergedTableRaw.txt -out $home/results/coverage/count-table-deeptools.merged.npz
