#!/bin/bash
#SBATCH --job-name=CAGE.bam2bed
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=5
#SBATCH --export=all
#SBATCH -p long



srx=("SRX882919" "SRX882920" "SRX882921" "SRX882922" "SRX882923" "SRX882924" "SRX882925" "SRX882926" "SRX882927" "SRX882928" "SRX882929" "SRX882930" "SRX882931" "SRX882932" "SRX882933" "SRX882934" "SRX882935" "SRX882936" "SRX882937" "SRX882938" "SRX882939" "SRX882940" "SRX882941")


for i in ${srx[@]}
do
  echo $i
  bedtools bamtobed -i $i".bam" > $i".bed"
done
