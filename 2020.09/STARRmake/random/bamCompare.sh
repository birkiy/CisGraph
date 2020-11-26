#!/bin/bash
#SBATCH --job-name=starr.compare
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task 20
#SBATCH --export=all
#SBATCH -p long
home=/groups/lackgrp/ll_members/berkay/STARRbegin/results

bamCompare \
  -b1 $home/mapping/processed/GR.0h.dex.merged.final.bam \
  -b2 $home/mapping/processed/input.GSE114063.merged.final.bam \
  -o $home/bigwig/GR.0h.dex.plasmid.merged.bigWig \
  --extendReads -p 20

bamCompare \
  -b1 $home/mapping/processed/GR.1h.dex.merged.final.bam \
  -b2 $home/mapping/processed/input.GSE114063.merged.final.bam \
  -o $home/bigwig/GR.1h.dex.plasmid.merged.bigWig \
  --extendReads -p 20

bamCompare \
  -b1 $home/mapping/processed/GR.4h.dex.merged.final.bam \
  -b2 $home/mapping/processed/input.GSE114063.merged.final.bam \
  -o $home/bigwig/GR.4h.dex.plasmid.merged.bigWig \
  --extendReads -p 20

bamCompare \
  -b1 $home/mapping/processed/GR.8h.dex.merged.final.bam \
  -b2 $home/mapping/processed/input.GSE114063.merged.final.bam \
  -o $home/bigwig/GR.8h.dex.plasmid.merged.bigWig \
  --extendReads -p 20

bamCompare \
  -b1 $home/mapping/processed/GR.12h.dex.merged.final.bam \
  -b2 $home/mapping/processed/input.GSE114063.merged.final.bam \
  -o $home/bigwig/GR.12h.dex.plasmid.merged.bigWig \
  --extendReads -p 20
