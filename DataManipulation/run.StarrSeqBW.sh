#!/bin/bash
#SBATCH --job-name=vers2.bw
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=15
#SBATCH --export=all
#SBATCH -p long




mapping=/groups/lackgrp/ll_members/tunc/phd/ana-starrseq-lncap-lacklab/analysis/mapping

bigwig=/groups/lackgrp/ll_members/tunc/phd/ana-starrseq-lncap-lacklab/analysis/bigwig


eth=lncap-etoh-merged.bam
dht=lncap-dht-merged.bam

outData=/home/ualtintas/github/Data/CisGraph/Vers2.0

--pseudocount 0.001 \
--operation log2 \


bamCoverage --bam $mapping/$eth -o $outData/Features/StarrSeq/starrseq.etoh.BPM.bigWig --samFlagInclude 64 --extendReads --normalizeUsing BPM -p 15
bamCoverage --bam $mapping/$dht -o $outData/Features/StarrSeq/starrseq.dht.BPM.bigWig --samFlagInclude 64 --extendReads --normalizeUsing BPM -p 15
