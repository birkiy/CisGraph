#!/bin/bash
#SBATCH --job-name=star.unique
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=30G
#SBATCH --cpus-per-task=20
#SBATCH --export=all
#SBATCH -p long


# samtools merge -@ 20 -u input.pool.bam \
#    results/mapping/input.GSE114063.pool1.bam \
#    results/mapping/input.GSE114063.pool2.bam \
#    results/mapping/input.GSE114063.pool3.bam \
#    results/mapping/input.GSE114063.pool4.bam \
#    results/mapping/input.GSE114063.pool5.bam \
#    results/mapping/input.GSE114063.pool6.bam \
#    results/mapping/input.GSE114063.pool7.bam \
#    results/mapping/input.GSE114063.pool8.bam \
#    results/mapping/input.GSE114063.pool9.bam \
#    results/mapping/input.GSE114063.pool10.bam \
#    results/mapping/input.GSE114063.pool11.bam \
#    results/mapping/input.GSE114063.pool12.bam



samtools view -Sb input.pool.bam | \
samtools sort -n -@ 20 -m 15G | \
samtools fixmate -@ 20 -r -m - - | \
samtools sort -@ 20 -m 15G | \
samtools markdup -@ 20 - input.pool.final.bam

samtools view -f2 -F2048 -b input.pool.final.bam | \
bedtools bamtobed -bedpe > input.pool.bedpe


cut -f 1,2,6,7 input.pool.bedpe > input.pool.bed











input.pool.bam

sort -k1,1 -k2,2n --parallel=10

cat \

   results/coverage/unique.library.pool1.bedpe \
   results/coverage/unique.library.pool2.bedpe \
   results/coverage/unique.library.pool3.bedpe \
   results/coverage/unique.library.pool4.bedpe \
   results/coverage/unique.library.pool5.bedpe \
   results/coverage/unique.library.pool6.bedpe \
   results/coverage/unique.library.pool7.bedpe \
   results/coverage/unique.library.pool8.bedpe \
   results/coverage/unique.library.pool9.bedpe \
   results/coverage/unique.library.pool10.bedpe \
   results/coverage/unique.library.pool11.bedpe \
   results/coverage/unique.library.pool12.bedpe > results/coverage/unique.library.bed

sort -k1,1 -k2,2n --parallel=10 results/coverage/unique.library.bed > results/coverage/unique.library.sort.bed


| \
   bedtools sort -i - | bedtools merge -i - > \
   results/coverage/unique.library.sort.ne.bedpe
