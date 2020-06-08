#!/bin/bash
#SBATCH --job-name=vers2.bw
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=15
#SBATCH --export=all
#SBATCH -p long


idx=/home/ualtintas/genomeAnnotations/bowtieIdx/hg19




outPath=/home/ualtintas/LncapExo/CAGE



srx=(
"SRX882919"
"SRX882920"
"SRX882921"
"SRX882922"
"SRX882923"
"SRX882924"
"SRX882925"
"SRX882926"
"SRX882927"
"SRX882928"
"SRX882929"
"SRX882930"
"SRX882931"
"SRX882932"
)


for SRX in ${srx[@]}
do
  echo $SRX;
  fastq=$SRX/*.fastq.gz
  bowtie -p 12 -n 2 -k 1 -m 1 -f -S --best $idx -q $fastq | \
  samtools view - -@ 15 -F 4 -Sb | \
  samtools sort - -@ 15 -n | \
  samtools fixmate - - -m | \
  samtools sort - -@ 15 | \
  samtools markdup - - -r | \
  samtools view - -@ 15 -b -q 30 -o $outPath/$SRX".processed.bam";
  samtools index $outPath/$SRX".processed.bam";
done


samtools merge -@ 10 $outPath/cage.etoh.bam $outPath/SRX882919.processed.bam $outPath/SRX882926.processed.bam
samtools merge -@ 10 $outPath/cage.dht.bam $outPath/SRX882925.processed.bam $outPath/SRX882932.processed.bam
