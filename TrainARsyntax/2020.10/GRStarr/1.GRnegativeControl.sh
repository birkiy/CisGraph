


home=/groups/lackgrp/ll_members/berkay/STARRbegin/results/peaks

bedtools intersect -a $home/GSM2095218_DEX_3hr_GR_ChIPseq.Rep1_peaks.bed -b $home/GSM2095219_DEX_3hr_GR_ChIPseq.Rep2_peaks.bed -u > $home/common.GR.peaks.bed


# Results for GCR_HUMAN.H11MO.0.A scan against hg19: 83712 hits
https://ccg.epfl.ch/pwmtools/wwwtmp/pwmscan_hg19_31609_43612.bed
bedtools intersect -a $home/pwmscan_hg19_31609_43612.bed -b $home/common.GR.peaks.bed -v > $home/negativeControlGR.bed

shuf -n 5000 $home/negativeControlGR.bed > $home/negativeControlGR.random.bed

cut -f 1,2,3 $home/negativeControlGR.random.bed | awk -v s=350 '{print $1, $2-s, $3+s, "negativeControlGR."NR}' | tr ' ' '\t' > $home/negativeControlGR.final.bed



bedtools sort -i $home/significantNegativeControlGR.bed | bedtools merge -i - > $home/significantNegativeControlGR.merged.bed
bedtools intersect -a $home/significantNegativeControlGR.merged.bed -b  $home/GSM2095218_DEX_3hr_GR_ChIPseq.Rep1_peaks.bed $home/GSM2095219_DEX_3hr_GR_ChIPseq.Rep2_peaks.bed -u






home=/groups/lackgrp/ll_members/berkay/STARRbegin/results/peaks


cat $home/GSM2095218_DEX_3hr_GR_ChIPseq.Rep1_peaks.bed $home/GSM2095219_DEX_3hr_GR_ChIPseq.Rep2_peaks.bed | \
  bedtools sort -i - | bedtools merge -i - | bedtools intersect -a $home/pwmscan_hg19_31609_43612.bed -b - -v > $home/negativeControlGR.bed


shuf -n 5000 $home/negativeControlGR.bed > $home/negativeControlGR.random.bed

cut -f 1,2,3 $home/negativeControlGR.random.bed | awk -v s=350 '{print $1, $2-s, $3+s, "negativeControlGR."NR}' | tr ' ' '\t' > $home/negativeControlGR.final.bed

cat $home/common.GR.peaks.bed $home/negativeControlGR.final.bed > $home/totalGRE.bed
