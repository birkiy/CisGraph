#!/bin/bash
#SBATCH --job-name=gro.vers2
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=3G
#SBATCH --cpus-per-task=12
#SBATCH --export=all
#SBATCH -p long


mapping=/groups/lackgrp/ll_members/common/LNCaP_bigwigs/groseq_analysis/mapping
outData=/home/ualtintas/github/Data/CisGraph/Vers2.0
bamCoverage -b $mapping/groseq_dht.bam --filterRNAstrand forward --scaleFactor 1 -o $outData/Features/GroSeq/groseq.dht.BPM.minus.bigWig  --normalizeUsing BPM -p 12
bamCoverage -b $mapping/groseq_dht.bam --filterRNAstrand reverse --scaleFactor 1 -o $outData/Features/GroSeq/groseq.dht.BPM.plus.bigWig --normalizeUsing BPM -p 12


bamCoverage -b $mapping/groseq_dmso.ba --filterRNAstrand forward --scaleFactor 1 -o $outData/Features/GroSeq/groseq.dmso.minus.bigWig  --normalizeUsing BPM -p 12
bamCoverage -b $mapping/groseq_dmso.ba --filterRNAstrand reverse --scaleFactor 1 -o $outData/Features/GroSeq/groseq.dmso.plus.bigWig --normalizeUsing BPM -p 12
