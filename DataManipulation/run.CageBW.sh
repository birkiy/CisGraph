#!/bin/bash
#SBATCH --job-name=vers2.bw
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=15
#SBATCH --export=all
#SBATCH -p long




mapping=/home/ualtintas/LncapExo/CAGE

bigwig=/groups/lackgrp/ll_members/tunc/phd/ana-starrseq-lncap-lacklab/analysis/bigwig


eth=cage.etoh.bam
dht=cage.dht.bam

outData=/home/ualtintas/github/Data/CisGraph/Vers2.0


bamCoverage --bam $mapping/$eth -o $outData/Features/Cage/cage.+.etoh.BPM.bigWig --normalizeUsing BPM -p 15 --filterRNAstrand forward
bamCoverage --bam $mapping/$eth -o $outData/Features/Cage/cage.-.etoh.BPM.bigWig --normalizeUsing BPM -p 15 --filterRNAstrand reverse



bamCoverage --bam $mapping/$dht -o $outData/Features/Cage/cage.+.dht.BPM.bigWig --normalizeUsing BPM -p 15 --filterRNAstrand forward
bamCoverage --bam $mapping/$dht -o $outData/Features/Cage/cage.-.dht.BPM.bigWig --normalizeUsing BPM -p 15 --filterRNAstrand reverse
