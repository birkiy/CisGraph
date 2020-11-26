

cat GSM803357_hg19_wgEncodeHaibTfbsA549GrPcr1xDex500pmPkRep1.broadPeak \
  GSM803357_hg19_wgEncodeHaibTfbsA549GrPcr1xDex500pmPkRep2.broadPeak \
  GSM803358_hg19_wgEncodeHaibTfbsA549GrPcr1xDex50nmPkRep1.broadPeak \
  GSM803358_hg19_wgEncodeHaibTfbsA549GrPcr1xDex50nmPkRep2.broadPeak \
  GSM803371_hg19_wgEncodeHaibTfbsA549GrPcr2xDex100nmPkRep1.broadPeak \
  GSM803371_hg19_wgEncodeHaibTfbsA549GrPcr2xDex100nmPkRep2.broadPeak \
  | bedtools sort | bedtools merge > GR.bed

bedtools intersect -a pwmscan_hg19_39036_18370.bed -b GR.bed > negativeControlGR.bed
