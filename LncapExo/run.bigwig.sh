

bamCoverage -b SRX882911.bam --scaleFactor 1 -o miRNA.0h.DHT.bigWig  --normalizeUsing BPM -p 12
bamCoverage -b SRX882912.bam --scaleFactor 1 -o miRNA.6h.DHT.bigWig  --normalizeUsing BPM -p 12
bamCoverage -b SRX882913.bam --scaleFactor 1 -o miRNA.24h.DHT.bigWig  --normalizeUsing BPM -p 12



#!/bin/bash
#SBATCH --job-name=CAGE.bw
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=15
#SBATCH --export=all
#SBATCH -p long

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
name=("LNCaP.0h.rep1"
"LNCaP.3h.rep1"
"LNCaP.6h.rep1"
"LNCaP.12h.rep1"
"LNCaP.18h.rep1"
"LNCaP.48h.rep1"
"LNCaP.72h.rep1"
"LNCaP.0h.rep2"
"LNCaP.3h.rep2"
"LNCaP.6h.rep2"
"LNCaP.12h.rep2"
"LNCaP.18h.rep2"
"LNCaP.48h.rep2"
"LNCaP.72h.rep2"
"LNCaP.0h"
"LNCaP.Ethanol.6h"
"LNCaP.DHT.6h"
"LNCaP.DHT.+.Bicaltamide.6h"
"LNCaP.Bicaltamide.6h"
"LNCaP.Ethanol.24h"
"LNCaP.DHT.24h"
"LNCaP.DHT.+.Bicaltamide.24h"
"LNCaP.Bicaltamide.24h"
)




for ((i = 0; i < 23; i++))
do
  bamCoverage -b ${srx[i]}".bam" --scaleFactor 1 -o ${name[i]}".bigWig"  --normalizeUsing BPM -p 12
done









tsP=~/genomeAnnotations/Regions/TSS.hg19.+.bed
tsM=~/genomeAnnotations/Regions/TSS.hg19.-.bed


computeMatrix reference-point -S CAGE/*.bigWig -R $tsP $tsM --referencePoint center -a 3000 -b 3000 -o mat/CAGE.TSS.3.mat.gz -p 30

plotHeatmap -m mat/CAGE.TSS.3.mat.gz -out CAGE.TSS.3.pdf --heatmapWidth 10 --heatmapHeight 50 --colorMap "YlGnBu"




con=~/ARBSs/regions/cons-arbs.bed
ind=~/ARBSs/regions/ind-arbs.bed
non=~/ARBSs/regions/Non-Active-ARBS.bed



computeMatrix reference-point -S CAGE/*.bigWig -R $con $ind $non --referencePoint center -a 3000 -b 3000 -o mat/CAGE.ARBS.3.mat.gz -p 30

plotHeatmap -m mat/CAGE.ARBS.3.mat.gz -out CAGE.ARBS.3.pdf --heatmapWidth 10 --heatmapHeight 50 --colorMap "YlGnBu"
