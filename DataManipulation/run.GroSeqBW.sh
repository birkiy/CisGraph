#!/bin/bash
#SBATCH --job-name=gro.vers2
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=3G
#SBATCH --cpus-per-task=12
#SBATCH --export=all
#SBATCH -p long


mapping=/groups/lackgrp/ll_members/common/LNCaP_bigwigs/groseq_analysis/mapping
outData=/home/ualtintas/github/Data/CisGraph/Vers2.0
bl=$outData/Blacklist/ENCFF001TDO.bed
blO=$outData/Blacklist/ENCFF001TDO.merged.bed
bedtools merge -i $bl -o collapse -c 4 -delim "-" > $blO

bamCoverage -b $mapping/groseq_dht.bam --filterRNAstrand forward -o $outData/Features/GroSeq/groseq.dht.-.BPM.bigWig  --normalizeUsing BPM -p 12 -bl $blO
bamCoverage -b $mapping/groseq_dht.bam --filterRNAstrand reverse -o $outData/Features/GroSeq/groseq.dht.+.BPM.bigWig --normalizeUsing BPM -p 12 -bl $blO


bamCoverage -b $mapping/groseq_dmso.bam --filterRNAstrand forward -o $outData/Features/GroSeq/groseq.dmso.-.BPM.bigWig  --normalizeUsing BPM -p 12 -bl $blO
bamCoverage -b $mapping/groseq_dmso.bam --filterRNAstrand reverse -o $outData/Features/GroSeq/groseq.dmso.+.BPM.bigWig --normalizeUsing BPM -p 12 -bl $blO
