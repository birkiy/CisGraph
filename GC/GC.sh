bedtools getfasta -fi ~/genomeAnnotations/hg19/hg19.fa -bed ../promoters_ann_5kb.bed -fo promoters.fasta

cut -f1,2,3,4 GSE86832_bigTable.tsv > GCmet.tsv

cut -f2,3,4 GCmet.tsv > end
cut -f1,2 GCmet.tsv > start

paste start end > GCmet.bed

awk -F'\t' '{ if(($5 != 0 ) && ($4 != 0)) {printf "%s\t%s\t%s\t%s.%s\t%s\t%s\n", $1, $2, $3, $1, NR, $4, $5}}' GCmet.bed > GC.bed


bedtools intersect -a GC.bed -b cons-arbs.bed ind-arbs.bed Non-Active-ARBS.bed > GCint.bed
