



bedtools intersect -a AR.dht.rep1.sorted.peaks.narrowPeak -b AR.dht.rep2.sorted.peaks.narrowPeak



head AR.dht.rep1.sorted.peaks.narrowPeak -n 5000 > AR.dht.rep1.head5k.narrowPeak
head AR.dht.rep2.sorted.peaks.narrowPeak -n 5000 > AR.dht.rep2.head5k.narrowPeak

head FOXA1.dht.rep1.sorted.peaks.narrowPeak -n 5000 > FOXA1.dht.rep1.head5k.narrowPeak
head FOXA1.dht.rep2.sorted.peaks.narrowPeak -n 5000 > FOXA1.dht.rep2.head5k.narrowPeak


head HOXB13.dht.rep1.sorted.peaks.narrowPeak -n 5000 > HOXB13.dht.rep1.head5k.narrowPeak




#!/bin/bash
#SBATCH --job-name=Chip.download
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=30
#SBATCH --export=all
#SBATCH -p long


bigwigCompare \
-b1 /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/bigwig/AR.dht.rep1.+.5end.bigWig \
-b2 /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/bigwig/AR.dht.rep2.+.5end.bigWig \
-o /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/bigwig/AR.dht.+.5end.bigWig \
--operation=add -bl /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed

bigwigCompare \
-b1 /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/bigwig/FOXA1.dht.rep1.+.5end.bigWig \
-b2 /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/bigwig/FOXA1.dht.rep2.+.5end.bigWig \
-o /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/bigwig/FOXA1.dht.+.5end.bigWig \
--operation=add -bl /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed



bigwigCompare \
-b1 /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/bigwig/AR.dht.rep1.-.5end.bigWig \
-b2 /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/bigwig/AR.dht.rep2.-.5end.bigWig \
-o /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/bigwig/AR.dht.+.5end.bigWig \
--operation=add -bl /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed

bigwigCompare \
-b1 /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/bigwig/FOXA1.dht.rep1.-.5end.bigWig \
-b2 /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/bigwig/FOXA1.dht.rep2.-.5end.bigWig \
-o /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/bigwig/FOXA1.dht.-.5end.bigWig \
--operation=add -bl /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed




bigwigCompare \
-b1 /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/bigwig/AR.dht.rep1.+.5end.bigWig \
-b2 /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/bigwig/AR.dht.rep2.+.5end.bigWig \
-o /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/bigwig/AR.dht.+.5end.bigWig \
--operation=add -bl /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed

bigwigCompare \
-b1 /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/bigwig/FOXA1.dht.rep1.+.5end.bigWig \
-b2 /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/bigwig/FOXA1.dht.rep2.+.5end.bigWig \
-o /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/bigwig/FOXA1.dht.+.5end.bigWig \
--operation=add -bl /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed









bigwigCompare \
-b1 /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/bigwig/control.GSE83860.rep1.-.5end.bigWig \
-b2 /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/bigwig/control.GSE83860.rep2.-.5end.bigWig \
-o /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/bigwig/control.GSE83860.-.5end.bigWig \
--operation=add -bl /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed

bigwigCompare \
-b1 /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/bigwig/control.GSE83860.rep1.-.5end.bigWig \
-b2 /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/bigwig/control.GSE83860.rep2.-.5end.bigWig \
-o /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/bigwig/control.GSE83860.-.5end.bigWig \
--operation=add -bl /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed




bigwigCompare \
-b1 /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/bigwig/control.GSE83860.rep1.+.5end.bigWig \
-b2 /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/bigwig/control.GSE83860.rep2.+.5end.bigWig \
-o /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/bigwig/control.GSE83860.+.5end.bigWig \
--operation=add -bl /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed

bigwigCompare \
-b1 /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/bigwig/control.GSE83860.rep1.+.5end.bigWig \
-b2 /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/bigwig/control.GSE83860.rep2.+.5end.bigWig \
-o /groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/bigwig/control.GSE83860.+.5end.bigWig \
--operation=add -bl /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed
