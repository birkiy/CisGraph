# exclude a large portion of the training set
exclude_chr=["chrX","chrY","chr5","chr6","chr7","chr10","chr14","chr11","chr13","chr12","chr15"]
valid_chr = ['chr2']
test_chr = ['chr1', 'chr8', 'chr9',
            'chr3', 'chr4']
seq_width = 200
n_dil_layers = 3
bpnet_data.interval_augmentation_shift = 100
train.seed = 42








# setup all the file paths

exp_dir="/groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/train"
model_dir=$exp_dir/'output'
contrib_file=$model_dir/'contrib.deeplift.h5'
contrib_null_file=$model_dir/'contrib.deeplift.null.h5'
modisco_dir=$model_dir/'modisco'


bpnet train dataspec.yml --premade=bpnet9 --config=config.gin . --override='train.epochs=10' --in-memory --num-workers 16 --run-id='exampleRun'



#!/bin/bash
#SBATCH --job-name=bpnet.bw
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task 30
#SBATCH --export=all
#SBATCH -p long


bamCoverage --bam /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/mapping/processed/AR.dht.rep1.final.bam \
            -of bigwig -o /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/bigwig/AR.dht.rep1.forward.ex.bigWig \
            --filterRNAstrand forward \
            --extendReads 150 \
            --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed \
            -p 30 --binSize 1

bamCoverage --bam /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/mapping/processed/AR.dht.rep1.final.bam \
            -of bigwig -o /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/bigwig/AR.dht.rep1.reverse.ex.bigWig \
            --filterRNAstrand reverse \
            --extendReads 150 \
            --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed \
            -p 30 --binSize 1



bamCoverage --bam /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/mapping/processed/FOXA1.dht.rep1.final.bam \
            -of bigwig -o /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/bigwig/FOXA1.dht.rep1.forward.ex.bigWig \
            --filterRNAstrand forward \
            --extendReads 150 \
            --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed \
            -p 30 --binSize 1

bamCoverage --bam /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/mapping/processed/FOXA1.dht.rep1.final.bam \
            -of bigwig -o /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/bigwig/FOXA1.dht.rep1.reverse.ex.bigWig \
            --filterRNAstrand reverse \
            --extendReads 150 \
            --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed \
            -p 30 --binSize 1


bamCoverage --bam /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/mapping/processed/control.GSE83860.rep1.final.bam \
            -of bigwig -o /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/bigwig/control.GSE83860.rep1.forward.ex.bigWig \
            --filterRNAstrand forward \
            --extendReads 150 \
            --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed \
            -p 30 --binSize 1

bamCoverage --bam /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/mapping/processed/control.GSE83860.rep1.final.bam \
            -of bigwig -o /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/bigwig/control.GSE83860.rep1.reverse.ex.bigWig \
            --filterRNAstrand reverse \
            --extendReads 150 \
            --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed \
            -p 30 --binSize 1
