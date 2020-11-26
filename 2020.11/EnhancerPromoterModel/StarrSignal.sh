#!/bin/bash
#SBATCH --job-name=starrSignal
#SBATCH --time=23:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=20
#SBATCH --export=all
#SBATCH -p long


bam=/groups/lackgrp/ll_members/tunc/phd/ana-starrseq-lncap-lacklab/analysis/mapping


home=/groups/lackgrp/ll_members/berkay/enhancerPromoterModel/StarrSignal

multiBamSummary BED-file \
  --BED $home/ARBS.bed \
  --bamfiles \
    $bam/lncap-dht-merged.bam \
    $bam/lncap-etoh-merged.bam \
    $bam/starrseq-lncap-lib3-cleaned.bam \
  --extendReads --samFlagInclude 64 \
  --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed \
  --centerReads \
  -p 20 \
  --outRawCounts AR.starr.signal.tsv -out AR.starr.signal.npz


tss=~/genomeAnnotations/Regions/TSS.hg19.Idx.bed

con=~/ARBSs/regions/cons-arbs.bed
ind=~/ARBSs/regions/ind-arbs.bed
non=~/ARBSs/regions/Non-Active-ARBS.bed
nAR=~/ARBSs/regions/negativeControl.ARBS.bed


cat $con $ind $non $nAR > ARBS.bed
