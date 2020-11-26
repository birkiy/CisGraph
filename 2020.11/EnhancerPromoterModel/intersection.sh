

home=/groups/lackgrp/ll_members/berkay/enhancerPromoterModel/intersect

sureCounts=/groups/lackgrp/ll_members/berkay/enhancerPromoterModel/supplementary_dataset1/SuRE-peaks_K562.45.55_raw_sep_globalLambda.annotated_LP160616.txt.gz



paste <(zcat $sureCounts | cut -f2,3,4)  <(zcat $sureCounts | cut -f1,6,7) | tail -n+2 > $home/SuRE.bed



bedtools intersect -wa -wb -a $home/SuRE.bed -b ~/ARBSs/regions/cons-arbs.bed > $home/con.sure.bed

bedtools intersect -wa -wb -a $home/SuRE.bed -b ~/ARBSs/regions/ind-arbs.bed > $home/ind.sure.bed

bedtools intersect -wa -wb -a $home/SuRE.bed -b ~/ARBSs/regions/Non-Active-ARBS.bed > $home/non.sure.bed

bedtools intersect -wa -wb -a $home/SuRE.bed -b ~/ARBSs/regions/negativeControl.ARBS.bed > $home/nAR.sure.bed

bedtools intersect -wa -wb -a $home/SuRE.bed -b ~/genomeAnnotations/Regions/TSS.hg19.Idx.bed > $home/pro.sure.bed
