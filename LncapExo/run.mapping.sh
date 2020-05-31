#!/bin/bash
#SBATCH --job-name=miR.map
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=15
#SBATCH --export=all
#SBATCH -p long


# bowtie2-build -f hg19.fa hg19
# bowtie-build hg19.fa hg19


home=/home/ualtintas/LncapExo/miRNA/
srx=("SRX882911"  "SRX882912" "SRX882913")

idx=/home/ualtintas/genomeAnnotations/bowtieIdx/hg19

for i in ${srx[@]}
do
        filename=`basename $i`
        fastq=$filename"/*.fastq.gz"
        echo $filename
        # bowtie --wrapper basic-0 --wrapper basic-0 -p 12 -v 1 -k 1 -m 1 -f -S --best -e 99999 $idx -q $fastq $filename".sam"
        bowtie -p 12 -v 1 -k 1 -m 1 -f -S --best ~/tools/bowtie_index/hg19 -q $i  | samtools view -Sb1 - | samtools sort -m 4G -@ 12 - -o $filename\.bam
        samtools index $filename\.bam
done



#!/bin/bash
#SBATCH --job-name=CAGE.map
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=15
#SBATCH --export=all
#SBATCH -p long


home=/home/ualtintas/LncapExo/CAGE/
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
"SRX882933"
"SRX882934"
"SRX882935"
"SRX882936"
"SRX882937"
"SRX882938"
"SRX882939"
"SRX882940"
"SRX882941"
)
idx=/home/ualtintas/genomeAnnotations/bowtieIdx/hg19

for i in ${srx[@]}
do
        filename=`basename $i`
        fastq=$filename"/*.fastq.gz"
        echo $filename
        # bowtie --wrapper basic-0 --wrapper basic-0 -p 12 -v 1 -k 1 -m 1 -f -S --best -e 99999 $idx -q $fastq $filename".sam"
        bowtie -p 12 -v 1 -k 1 -m 1 -f -S --best $idx -q $fastq  | samtools view -Sb1 - | samtools sort -m 4G -@ 12 - -o $filename\.bam
        samtools index $filename\.bam
done
