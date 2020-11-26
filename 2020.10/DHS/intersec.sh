


home=/groups/lackgrp/ll_members/berkay/DHS


tail -n +2 $home/rawData/DHS_Index_and_Vocabulary_hg19_WM20190703.txt | awk -F'\t' '{print $0"\t"NR}' - | \
  bedtools intersect -a - -b ~/ARBSs/regions/cons-arbs.bed > $home/creOverARBS/con.cre.bed

tail -n +2 $home/rawData/DHS_Index_and_Vocabulary_hg19_WM20190703.txt|  awk -F'\t' '{print $0"\t"NR}' - | \
  bedtools intersect -a - -b ~/ARBSs/regions/ind-arbs.bed > $home/creOverARBS/ind.cre.bed

tail -n +2 $home/rawData/DHS_Index_and_Vocabulary_hg19_WM20190703.txt|  awk -F'\t' '{print $0"\t"NR}' - | \
  bedtools intersect -a - -b ~/ARBSs/regions/Non-Active-ARBS.bed > $home/creOverARBS/non.cre.bed

tail -n +2 $home/rawData/DHS_Index_and_Vocabulary_hg19_WM20190703.txt|  awk -F'\t' '{print $0"\t"NR}' - | \
  bedtools intersect -a - -b ~/ARBSs/regions/negativeControl.ARBS.bed > $home/creOverARBS/nAR.cre.bed








sort -n <(cut -f7 $home/creOverARBS/con.cre.bed) \
  <(cut -f7 $home/creOverARBS/ind.cre.bed) \
  <(cut -f7 $home/creOverARBS/non.cre.bed) \
  <(cut -f7 $home/creOverARBS/nAR.cre.bed) | \
  uniq > $home/creOverARBS/ARBS.index.txt

awk 'NR==FNR{data[$1]; next}FNR in data' \
  $home/creOverARBS/ARBS.index.txt \
  $home/rawData/dat_bin_FDR01_hg19.txt \
  > $home/creOverARBS/indexed.binary.txt
