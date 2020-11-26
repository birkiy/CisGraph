#!/bin/bash
#SBATCH --job-name=starr.count
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task 64
#SBATCH --export=all
#SBATCH -p long


home=/groups/lackgrp/ll_members/berkay/StarrPipe/results/mapping/each/coorsorted

multiBamSummary BED-file \
  --BED /groups/lackgrp/ll_members/berkay/GenomicRegions/GR/common.GR.peaks.bed \
  /groups/lackgrp/ll_members/berkay/GenomicRegions/GR/negativeControlGR.final.bed \
  --bamfiles \
  $home/A549.GR.dex.0h.rep1.bam \
  $home/A549.GR.dex.0h.rep2.bam \
  $home/A549.GR.dex.0h.rep3.bam \
  $home/A549.GR.dex.0h.rep4.bam \
  $home/A549.GR.dex.0h.rep5.bam \
  $home/A549.GR.dex.12h.rep1.bam \
  $home/A549.GR.dex.12h.rep2.bam \
  $home/A549.GR.dex.12h.rep3.bam \
  $home/A549.GR.dex.12h.rep4.bam \
  $home/A549.GR.dex.12h.rep5.bam \
  $home/A549.GR.dex.1h.rep1.bam \
  $home/A549.GR.dex.1h.rep2.bam \
  $home/A549.GR.dex.1h.rep3.bam \
  $home/A549.GR.dex.1h.rep4.bam \
  $home/A549.GR.dex.1h.rep5.bam \
  $home/A549.GR.dex.4h.rep1.bam \
  $home/A549.GR.dex.4h.rep2.bam \
  $home/A549.GR.dex.4h.rep3.bam \
  $home/A549.GR.dex.4h.rep4.bam \
  $home/A549.GR.dex.4h.rep5.bam \
  $home/A549.GR.dex.8h.rep1.bam \
  $home/A549.GR.dex.8h.rep2.bam \
  $home/A549.GR.dex.8h.rep3.bam \
  $home/A549.GR.dex.8h.rep4.bam \
  $home/A549.GR.dex.8h.rep5.bam \
  $home/GSE114063.10.control.bam \
  $home/GSE114063.11.control.bam \
  $home/GSE114063.12.control.bam \
  $home/GSE114063.1.control.bam \
  $home/GSE114063.2.control.bam \
  $home/GSE114063.3.control.bam \
  $home/GSE114063.4.control.bam \
  $home/GSE114063.5.control.bam \
  $home/GSE114063.6.control.bam \
  $home/GSE114063.7.control.bam \
  $home/GSE114063.8.control.bam \
  $home/GSE114063.9.control.bam \
  --extendReads \
  --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed \
  -p 64 \
  --outRawCounts results/coverage/each/multiBamSum/multiBamSum.tsv  \
  -out results/coverage/each/multiBamSum/multiBamSum.npz



home=/groups/lackgrp/ll_members/berkay/StarrPipe/results/mapping/merged/namesorted

multiBamSummary BED-file \
  --BED /groups/lackgrp/ll_members/berkay/GenomicRegions/GR/common.GR.peaks.bed \
  /groups/lackgrp/ll_members/berkay/GenomicRegions/GR/negativeControlGR.final.bed \
  --bamfiles \
    $home/A549.GR.dex.0h.bam \
    $home/A549.GR.dex.1h.bam \
    $home/A549.GR.dex.4h.bam \
    $home/A549.GR.dex.8h.bam \
    $home/A549.GR.dex.12h.bam \
    $home/GSE114063.control.bam \
  --extendReads \
  --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed \
  --centerReads \
  -p 64 \
  --outRawCounts results/coverage/merged/multiBamSum/multiBamSum.tsv  \
  -out results/coverage/merged/multiBamSum/multiBamSum.npz





plotPCA --transpose --in results/coverage/each/multiBamSum/multiBamSum.npz \
  -o results/plots/multiBamSum.each.PCA.transposed.pdf

plotPCA --in results/coverage/each/multiBamSum/multiBamSum.npz \
  -o results/plots/multiBamSum.each.PCA.pdf --labels "0h" "1h" "4h" "8h" "12h"

plotPCA --rowCenter --in results/coverage/each/multiBamSum/multiBamSum.npz \
  -o results/plots/multiBamSum.each.PCA.rowCenter.pdf
