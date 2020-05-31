
fi=/kuacc/users/ualtintas20/genomeAnnotations/hg19/hg19.fa

bed=/kuacc/users/ualtintas20/regions/LNCaP_DHT_DHS.bed

bedtools getfasta -fi $fi -bed $bed -fo ~/regions/sequences/creDHS.fasta -name 4
