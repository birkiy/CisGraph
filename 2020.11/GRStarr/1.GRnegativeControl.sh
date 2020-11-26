


https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2095nnn/GSM2095218/suppl/GSM2095218_DEX_3hr_GR_ChIPseq.Rep1_peaks.bed.gz
https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2095nnn/GSM2095219/suppl/GSM2095219_DEX_3hr_GR_ChIPseq.Rep2_peaks.bed.gz

bedtools intersect -a GSM2095218_DEX_3hr_GR_ChIPseq.Rep1_peaks.bed -b GSM2095219_DEX_3hr_GR_ChIPseq.Rep2_peaks.bed -u > common.GR.peaks.bed


# pwmtools => pwmscan
# HOCOMOCO v11 HUMAN TF
# Results for GCR_HUMAN.H11MO.0.A scan against hg19: 83712 hits
# https://ccg.epfl.ch/pwmtools/wwwtmp/pwmscan_hg19_14347_14769.bed


# LAST VERSION


cat GSM2095218_DEX_3hr_GR_ChIPseq.Rep1_peaks.bed GSM2095219_DEX_3hr_GR_ChIPseq.Rep2_peaks.bed | \
bedtools sort -i - | bedtools merge -i - | bedtools intersect -a pwmscan_hg19_33372_34079.bed -b - -v > negativeControlGR.bed


shuf -n 5000 negativeControlGR.bed > negativeControlGR.random.bed

cut -f 1,2,3 negativeControlGR.random.bed | awk -v s=350 '{print $1, $2-s, $3+s, "negativeControlGR."NR}' | tr ' ' '\t' > negativeControlGR.final.bed

cat common.GR.peaks.bed negativeControlGR.final.bed > totalGRE.bed



################################################################################




#
# bedtools intersect -a pwmscan_hg19_33372_34079.bed -b common.GR.peaks.bed -v > negativeControlGR.bed
#
# shuf -n 5000 negativeControlGR.bed > negativeControlGR.random.bed
#
# cut -f 1,2,3 negativeControlGR.random.bed | awk -v s=350 '{print $1, $2-s, $3+s, "negativeControlGR."NR}' | tr ' ' '\t' > negativeControlGR.final.bed

#############################################
#
# bedtools sort -i $home/significantNegativeControlGR.bed | bedtools merge -i - > $home/significantNegativeControlGR.merged.bed
# bedtools intersect -a $home/significantNegativeControlGR.merged.bed -b  $home/GSM2095218_DEX_3hr_GR_ChIPseq.Rep1_peaks.bed $home/GSM2095219_DEX_3hr_GR_ChIPseq.Rep2_peaks.bed -u