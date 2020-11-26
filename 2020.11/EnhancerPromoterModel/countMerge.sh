#!/bin/bash
#SBATCH --job-name=sureSignal
#SBATCH --time=23:00:00
#SBATCH --mem-per-cpu=300G
#SBATCH --cpus-per-task=1
#SBATCH --export=all
#SBATCH -p long


home=/groups/lackgrp/ll_members/berkay/enhancerPromoterModel/SuREsnpCount

cut -f1,2,3,4,10,11,12,13,14,15 $home/GSE128325_SuRE42_1* | tail -n+2 >> $home/GSE128325_SuRE42_1.count.bed

cut -f1,2,3,4,10,11,12,13,14,15 $home/GSE128325_SuRE42_2* | tail -n+2 > $home/GSE128325_SuRE42_2.count.bed


cut -f1,2,3,4,10,11,12,13,14,15 $home/GSE128325_SuRE43_1* | tail -n+2 > $home/GSE128325_SuRE43_1.count.bed

cut -f1,2,3,4,10,11,12,13,14,15 $home/GSE128325_SuRE43_2* | tail -n+2 > $home/GSE128325_SuRE43_2.count.bed


cut -f1,2,3,4,10,11,12,13,14,15 $home/GSE128325_SuRE44_1* | tail -n+2 > $home/GSE128325_SuRE44_1.count.bed

cut -f1,2,3,4,10,11,12,13,14,15 $home/GSE128325_SuRE44_2* | tail -n+2 > $home/GSE128325_SuRE44_2.count.bed


cut -f1,2,3,4,10,11,12,13,14,15 $home/GSE128325_SuRE45_1* | tail -n+2 > $home/GSE128325_SuRE45_1.count.bed

cut -f1,2,3,4,10,11,12,13,14,15 $home/GSE128325_SuRE45_2* | tail -n+2 > $home/GSE128325_SuRE45_2.count.bed



home=/groups/lackgrp/ll_members/berkay/enhancerPromoterModel/SuREsnpCount

tss=~/genomeAnnotations/Regions/TSS.hg19.Idx.bed

con=~/ARBSs/regions/cons-arbs.bed
ind=~/ARBSs/regions/ind-arbs.bed
non=~/ARBSs/regions/Non-Active-ARBS.bed
nAR=~/ARBSs/regions/negativeControl.ARBS.bed


bedtools intersect -wa -wb -a $home/GSE128325_SuRE42_1.count.bed -b $con > $home/ARBS/GSE128325_SuRE42_1.count.con.bed
bedtools intersect -wa -wb -a $home/GSE128325_SuRE42_2.count.bed -b $con > $home/ARBS/GSE128325_SuRE42_2.count.con.bed
bedtools intersect -wa -wb -a $home/GSE128325_SuRE43_1.count.bed -b $con > $home/ARBS/GSE128325_SuRE43_1.count.con.bed
bedtools intersect -wa -wb -a $home/GSE128325_SuRE43_2.count.bed -b $con > $home/ARBS/GSE128325_SuRE43_2.count.con.bed
bedtools intersect -wa -wb -a $home/GSE128325_SuRE44_1.count.bed -b $con > $home/ARBS/GSE128325_SuRE44_1.count.con.bed
bedtools intersect -wa -wb -a $home/GSE128325_SuRE44_2.count.bed -b $con > $home/ARBS/GSE128325_SuRE44_2.count.con.bed
bedtools intersect -wa -wb -a $home/GSE128325_SuRE45_1.count.bed -b $con > $home/ARBS/GSE128325_SuRE45_1.count.con.bed
bedtools intersect -wa -wb -a $home/GSE128325_SuRE45_2.count.bed -b $con > $home/ARBS/GSE128325_SuRE45_2.count.con.bed

bedtools intersect -wa -wb -a $home/GSE128325_SuRE42_1.count.bed -b $ind > $home/ARBS/GSE128325_SuRE42_1.count.ind.bed
bedtools intersect -wa -wb -a $home/GSE128325_SuRE42_2.count.bed -b $ind > $home/ARBS/GSE128325_SuRE42_2.count.ind.bed
bedtools intersect -wa -wb -a $home/GSE128325_SuRE43_1.count.bed -b $ind > $home/ARBS/GSE128325_SuRE43_1.count.ind.bed
bedtools intersect -wa -wb -a $home/GSE128325_SuRE43_2.count.bed -b $ind > $home/ARBS/GSE128325_SuRE43_2.count.ind.bed
bedtools intersect -wa -wb -a $home/GSE128325_SuRE44_1.count.bed -b $ind > $home/ARBS/GSE128325_SuRE44_1.count.ind.bed
bedtools intersect -wa -wb -a $home/GSE128325_SuRE44_2.count.bed -b $ind > $home/ARBS/GSE128325_SuRE44_2.count.ind.bed
bedtools intersect -wa -wb -a $home/GSE128325_SuRE45_1.count.bed -b $ind > $home/ARBS/GSE128325_SuRE45_1.count.ind.bed
bedtools intersect -wa -wb -a $home/GSE128325_SuRE45_2.count.bed -b $ind > $home/ARBS/GSE128325_SuRE45_2.count.ind.bed

bedtools intersect -wa -wb -a $home/GSE128325_SuRE42_1.count.bed -b $non > $home/ARBS/GSE128325_SuRE42_1.count.non.bed
bedtools intersect -wa -wb -a $home/GSE128325_SuRE42_2.count.bed -b $non > $home/ARBS/GSE128325_SuRE42_2.count.non.bed
bedtools intersect -wa -wb -a $home/GSE128325_SuRE43_1.count.bed -b $non > $home/ARBS/GSE128325_SuRE43_1.count.non.bed
bedtools intersect -wa -wb -a $home/GSE128325_SuRE43_2.count.bed -b $non > $home/ARBS/GSE128325_SuRE43_2.count.non.bed
bedtools intersect -wa -wb -a $home/GSE128325_SuRE44_1.count.bed -b $non > $home/ARBS/GSE128325_SuRE44_1.count.non.bed
bedtools intersect -wa -wb -a $home/GSE128325_SuRE44_2.count.bed -b $non > $home/ARBS/GSE128325_SuRE44_2.count.non.bed
bedtools intersect -wa -wb -a $home/GSE128325_SuRE45_1.count.bed -b $non > $home/ARBS/GSE128325_SuRE45_1.count.non.bed
bedtools intersect -wa -wb -a $home/GSE128325_SuRE45_2.count.bed -b $non > $home/ARBS/GSE128325_SuRE45_2.count.non.bed

bedtools intersect -wa -wb -a $home/GSE128325_SuRE42_1.count.bed -b $nAR > $home/ARBS/GSE128325_SuRE42_1.count.nAR.bed
bedtools intersect -wa -wb -a $home/GSE128325_SuRE42_2.count.bed -b $nAR > $home/ARBS/GSE128325_SuRE42_2.count.nAR.bed
bedtools intersect -wa -wb -a $home/GSE128325_SuRE43_1.count.bed -b $nAR > $home/ARBS/GSE128325_SuRE43_1.count.nAR.bed
bedtools intersect -wa -wb -a $home/GSE128325_SuRE43_2.count.bed -b $nAR > $home/ARBS/GSE128325_SuRE43_2.count.nAR.bed
bedtools intersect -wa -wb -a $home/GSE128325_SuRE44_1.count.bed -b $nAR > $home/ARBS/GSE128325_SuRE44_1.count.nAR.bed
bedtools intersect -wa -wb -a $home/GSE128325_SuRE44_2.count.bed -b $nAR > $home/ARBS/GSE128325_SuRE44_2.count.nAR.bed
bedtools intersect -wa -wb -a $home/GSE128325_SuRE45_1.count.bed -b $nAR > $home/ARBS/GSE128325_SuRE45_1.count.nAR.bed
bedtools intersect -wa -wb -a $home/GSE128325_SuRE45_2.count.bed -b $nAR > $home/ARBS/GSE128325_SuRE45_2.count.nAR.bed

bedtools intersect -wa -wb -a $home/GSE128325_SuRE42_1.count.bed -b $tss > $home/ARBS/GSE128325_SuRE42_1.count.tss.bed
bedtools intersect -wa -wb -a $home/GSE128325_SuRE42_2.count.bed -b $tss > $home/ARBS/GSE128325_SuRE42_2.count.tss.bed
bedtools intersect -wa -wb -a $home/GSE128325_SuRE43_1.count.bed -b $tss > $home/ARBS/GSE128325_SuRE43_1.count.tss.bed
bedtools intersect -wa -wb -a $home/GSE128325_SuRE43_2.count.bed -b $tss > $home/ARBS/GSE128325_SuRE43_2.count.tss.bed
bedtools intersect -wa -wb -a $home/GSE128325_SuRE44_1.count.bed -b $tss > $home/ARBS/GSE128325_SuRE44_1.count.tss.bed
bedtools intersect -wa -wb -a $home/GSE128325_SuRE44_2.count.bed -b $tss > $home/ARBS/GSE128325_SuRE44_2.count.tss.bed
bedtools intersect -wa -wb -a $home/GSE128325_SuRE45_1.count.bed -b $tss > $home/ARBS/GSE128325_SuRE45_1.count.tss.bed
bedtools intersect -wa -wb -a $home/GSE128325_SuRE45_2.count.bed -b $tss > $home/ARBS/GSE128325_SuRE45_2.count.tss.bed



home=/groups/lackgrp/ll_members/berkay/enhancerPromoterModel

tss=~/genomeAnnotations/Regions/TSS.hg19.Idx.bed

con=~/ARBSs/regions/cons-arbs.bed
ind=~/ARBSs/regions/ind-arbs.bed
non=~/ARBSs/regions/Non-Active-ARBS.bed
nAR=~/ARBSs/regions/negativeControl.ARBS.bed

bedtools intersect -wa -wb -a <(tail -n+2 $home/GSE78709_SuRE-counts-K562_B45_B55_LP170105.txt) -b $tss > GSE78709_SuRE-counts.tss.bed
bedtools intersect -wa -wb -a <(tail -n+2 $home/GSE78709_SuRE-counts-K562_B45_B55_LP170105.txt) -b $con > GSE78709_SuRE-counts.con.bed
bedtools intersect -wa -wb -a <(tail -n+2 $home/GSE78709_SuRE-counts-K562_B45_B55_LP170105.txt) -b $ind > GSE78709_SuRE-counts.ind.bed
bedtools intersect -wa -wb -a <(tail -n+2 $home/GSE78709_SuRE-counts-K562_B45_B55_LP170105.txt) -b $non > GSE78709_SuRE-counts.non.bed
bedtools intersect -wa -wb -a <(tail -n+2 $home/GSE78709_SuRE-counts-K562_B45_B55_LP170105.txt) -b $nAR > GSE78709_SuRE-counts.nAR.bed







 <( )) >
bedtools sort -i $home/GSE128325_SuRE42_2.count.bed | bedtools merge -i -
bedtools sort -i $home/GSE128325_SuRE43_1.count.bed | bedtools merge -i -
bedtools sort -i $home/GSE128325_SuRE43_2.count.bed | bedtools merge -i -
bedtools sort -i $home/GSE128325_SuRE44_1.count.bed | bedtools merge -i -
bedtools sort -i $home/GSE128325_SuRE44_2.count.bed | bedtools merge -i -
bedtools sort -i $home/GSE128325_SuRE45_1.count.bed | bedtools merge -i -
bedtools sort -i $home/GSE128325_SuRE45_2.count.bed | bedtools merge -i -


tail -n+2 $home/GSE128325_SuRE42_1.count.bed | \
  bedtools sort -i - | \
  bedtools merge -i - | \
  >  $home/GSE128325_SuRE42_1.count.merged.bed




cut -f1,2,3,4,10,11,12,13,14,15 sample.txt | tail -n+2 | bedtools merge -i - | head
