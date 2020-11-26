#!/bin/bash
#SBATCH --job-name=star.snakemake
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=40
#SBATCH --export=all


snakemake --cores 40

head -n 10000 deneme.sam > denemeSub.sam
cat <(samtools view -H denemeSub.sam) <(samtools view -q 30 -F 1804 denemeSub.sam | awk '$6 !~ "I|D"' - ) | \
        samtools sort --threads 5 -m 10G -n - | \
        samtools fixmate -m --threads 5 - - | \
        samtools sort --threads 5 -m 10G - | \
        samtools markdup -r - - | \
        samtools view -b - > deneme.bam
      samtools index -@ 5 deneme.bam deneme.bam.bai



cat <(samtools view -H deneme2.sam) <(samtools view -q 30 -F 1804 deneme2.sam | awk '$6 !~ "I|D"' - ) | \
        samtools sort -@ 5 -m 10G -n - | \
        samtools fixmate -m -@ 5 - - | \
        samtools sort -@ 5 -m 10G - | \
        samtools markdup -r - - | \
        samtools view -b - > deneme2.bam


gzip -dc rawData/AR.dht.rep1.fastq.gz | bowtie --chunkmbs 512 -k 1 -m 1 -v 2 --best --strata "/home/ualtintas/genomeAnnotations/bowtieIdx/hg19" --threads 5 -q - -S deneme.sam


gzip -dc rawData/AR.dht.rep1.fastq.gz | bowtie --chunkmbs 512 -k 1 -m 1 -v 2 --best --strata "/home/ualtintas/genomeAnnotations/bowtieIdx/hg19" --threads 5 -q - -S deneme2.sam



bedtools genomecov -5 -bg -strand + -ibam deneme.bam | bedtools sort > deneme.bedGraph

bedtools genomecov -5 -bg -strand + -ibam deneme.bam | bedtools sort >deneme.bedGraph
bedGraphToBigWig deneme.bedGraph /home/ualtintas/genomeAnnotations/hg19.chrom.sizes deneme.bigWig

bedGraphToBigWig deneme.bedGraph /home/ualtintas/genomeAnnotations/hg19.chrom.sizes deneme.bigWig

snakemake --dag analysis/bigwig/AR.dht.rep1.reverse.5end.bigWig | dot -Tsvg > dag.svg

bedGraphToBigWig analysis/bigwig/TRIM28.r1881.rep2.forward.5end.bedGraph /home/ualtintas/genomeAnnotations/hg19.chrom.sizes analysis/bigwig/TRIM28.r1881.rep2.forward.5end.bigWig





bedtools genomecov -5 -bg -strand + -ibam analysis/mapping/WDHD1.veh.rep2.final.bam | sort -k1,1 -k2,2n > analysis/bigwig/WDHD1.veh.rep2.forward.5end.bedGraph
wc -l analysis/bigwig/WDHD1.veh.rep2.forward.5end.bedGraph
bedGraphToBigWig analysis/bigwig/WDHD1.veh.rep2.forward.5end.bedGraph /home/ualtintas/genomeAnnotations/hg19.chrom.sizes WDHD1.veh.rep2.forward.5end.bw

bedGraphToBigWig analysis/bigwig/WDHD1.veh.rep2.forward.5end.bedGraph /home/ualtintas/genomeAnnotations/hg19.chrom.sizes WDHD1.veh.rep2.forward.5end.bw
http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
