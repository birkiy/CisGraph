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
inp=starrseq-lncap-lib3-cleaned.bam


outData=/home/ualtintas/github/Data/CisGraph/Vers2.0


bamCoverage --bam $mapping/$eth -o $outData/Features/StarrSeq/starrseq.etoh.BPM.bigWig --samFlagInclude 64 --extendReads --normalizeUsing BPM -p 15
bamCoverage --bam $mapping/$dht -o $outData/Features/StarrSeq/starrseq.dht.BPM.bigWig --samFlagInclude 64 --extendReads --normalizeUsing BPM -p 15
bamCoverage --bam $mapping/$inp -o $outData/Features/StarrSeq/starrseq.inp3.BPM.bigWig --samFlagInclude 64 --extendReads --normalizeUsing BPM -p 15



bigwigCompare \
-b1 $outData/Features/StarrSeq/starrseq.etoh.BPM.bigWig \
-b2 $outData/Features/StarrSeq/starrseq.inp3.BPM.bigWig \
--pseudocount 0.001 \
--operation log2 \
-o $outData/Features/StarrSeq/starrseq.etoh.plasmid.BPM.bigWig \
-p 62


bigwigCompare \
-b1 $outData/Features/StarrSeq/starrseq.dht.BPM.bigWig \
-b2 $outData/Features/StarrSeq/starrseq.inp3.BPM.bigWig \
--pseudocount 0.001 \
--operation log2 \
-o $outData/Features/StarrSeq/starrseq.dht.plasmid.BPM.bigWig \
-p 62




bigwigCompare \
-b1 $outData/Features/StarrSeq/starrseq.inp3.BPM.bigWig \
-b2 $outData/Features/StarrSeq/starrseq.etoh.BPM.bigWig \
--pseudocount 0.001 \
--operation log2 \
-o $outData/Features/StarrSeq/starrseq.etoh.enhancer.BPM.bigWig \
-p 62


bigwigCompare \
-b1 $outData/Features/StarrSeq/starrseq.inp3.BPM.bigWig \
-b2 $outData/Features/StarrSeq/starrseq.dht.BPM.bigWig \
--pseudocount 0.001 \
--operation log2 \
-o $outData/Features/StarrSeq/starrseq.dht.enhancer.BPM.bigWig \
-p 62
