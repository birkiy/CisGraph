#!/bin/bash
#SBATCH --job-name=miR.Process
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=15
#SBATCH --export=all
#SBATCH -p long


srx=("SRX882911"  "SRX882912" "SRX882913")

home=/home/ualtintas/LncapExo/miRNA

for i in ${srx[@]}
do
        SRX=`basename $i`
        echo $SRX

      	samtools view -@ 15 -b -o $home"/"$SRX".bam" $home"/"$SRX".sam";

      	samtools view -@ 15 -F 4 -b -o $home"/"$SRX".mapped.bam" $home"/"$SRX".bam";

      	samtools sort -@ 15 -n -o $home"/"$SRX".mapped.collated.bam" $home"/"$SRX".mapped.bam" ;
      	samtools fixmate -m -O bam $home"/"$SRX".mapped.collated.bam" $home"/"$SRX".mapped.collated.fixed.bam" ;
      	samtools sort -@ 15 -o $home"/"$SRX".mapped.collated.fixed.possorted.bam" $home"/"$SRX".mapped.collated.fixed.bam"
      	samtools markdup -r $home"/"$SRX".mapped.collated.fixed.possorted.bam" $home"/"$SRX".mapped.collated.fixed.possorted.markdup.bam";

      	samtools view -@ 15 -b -q 20 -o $home"/"$SRX".mapped.collated.fixed.possorted.markdup.mapq.bam" $home"/"$SRX".mapped.collated.fixed.possorted.markdup.bam"

      	samtools index -b $home"/"$SRX".mapped.collated.fixed.possorted.markdup.mapq.bam"  $home"/"$SRX".final.bam.bai" ;

      	rm $home"/"$SRX".mapped.bam"
        rm $home"/"$SRX".mapped.collated.bam"
        rm $home"/"$SRX".mapped.collated.fixed.bam"
        rm $home"/"$SRX".mapped.collated.fixed.possorted.bam"
        rm $home"/"$SRX".mapped.collated.fixed.possorted.markdup.bam"

        mv $home"/"$SRX".mapped.collated.fixed.possorted.markdup.mapq.bam" $home"/"$SRX"final.bam"
done
