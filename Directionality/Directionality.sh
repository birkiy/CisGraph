
#!/bin/bash
#SBATCH --job-name=star.pbs
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=3G
#SBATCH --cpus-per-task=12
#SBATCH --export=all
#SBATCH -p express

home=/groups/lackgrp/ll_members/common/LNCaP_bigwigs/groseq_analysis/mapping/

bamCoverage -b $home"groseq_dht.bam" --filterRNAstrand forward --scaleFactor 1 -o groseq_dht.BPM.minus.bigWig  --normalizeUsing BPM -p 12
bamCoverage -b $home"groseq_dht.bam" --filterRNAstrand reverse --scaleFactor 1 -o groseq_dht.BPM.plus.bigWig --normalizeUsing BPM -p 12


bamCoverage -b $home"groseq_dmso.bam" --filterRNAstrand forward --scaleFactor 1 -o groseq_dmso.BPM.minus.bigWig  --normalizeUsing BPM -p 12
bamCoverage -b $home"groseq_dmso.bam" --filterRNAstrand reverse --scaleFactor 1 -o groseq_dmso.BPM.plus.bigWig --normalizeUsing BPM -p 12




# RPKM
# minDht=/groups/lackgrp/ll_members/common/LNCaP_bigwigs/groseq_analysis/mapping/groseq_dht.minus.bw
# pluDht=/groups/lackgrp/ll_members/common/LNCaP_bigwigs/groseq_analysis/mapping/groseq_dht.plus.bw
#
#
# minDms=/groups/lackgrp/ll_members/common/LNCaP_bigwigs/groseq_analysis/mapping/groseq_dmso.minus.bw
# pluDms=/groups/lackgrp/ll_members/common/LNCaP_bigwigs/groseq_analysis/mapping/groseq_dmso.plus.bw


# BPM
minDht=~/BigWig/groseq_dht.BPM.minus.bigWig
pluDht=~/BigWig/groseq_dht.BPM.plus.bigWig

minDms=~/BigWig/groseq_dmso.BPM.minus.bigWig
pluDms=~/BigWig/groseq_dmso.BPM.plus.bigWig


# Promoters
tss=~/genomeAnnotations/Regions/TSS.hg19.Idx.bed
tsP=~/genomeAnnotations/Regions/TSS.hg19.+.bed
tsM=~/genomeAnnotations/Regions/TSS.hg19.-.bed

# ARBSs
con=~/ARBSs/regions/cons-arbs.bed
ind=~/ARBSs/regions/ind-arbs.bed
non=~/ARBSs/regions/Non-Active-ARBS.bed


coP=~/ARBSs/DirectionARBS/con.+.bed
inP=~/ARBSs/DirectionARBS/ind.+.bed
noP=~/ARBSs/DirectionARBS/non.+.bed

coM=~/ARBSs/DirectionARBS/con.-.bed
inM=~/ARBSs/DirectionARBS/ind.-.bed
noM=~/ARBSs/DirectionARBS/non.-.bed



# FirstPart
awk -F'\t' '{print $1"\t"($2 + 150) " ""\t"($3 - 500) " ""\t"$4"\t"$5
}' $tss > Parts/tss.FirstPart.bed

awk -F'\t' '{print $1"\t"($2 + 150) " ""\t"($3 - 500) " ""\t"$4"\t"$5
}' $tsP > Parts/tsP.FirstPart.bed

awk -F'\t' '{print $1"\t"($2 + 150) " ""\t"($3 - 500) " ""\t"$4"\t"$5
}' $tsM > Parts/tsM.FirstPart.bed

awk -F'\t' '{print $1"\t"$2"\t"($3 - 350) " ""\t"$4"\t"$5
}' $con > Parts/con.FirstPart.bed

awk -F'\t' '{print $1"\t"$2"\t"($3 - 350) " ""\t"$4"\t"$5
}' $ind > Parts/ind.FirstPart.bed

awk -F'\t' '{print $1"\t"$2"\t"($3 - 350) " ""\t"$4"\t"$5
}' $non > Parts/non.FirstPart.bed


# SecondPart
awk -F'\t' '{print $1"\t"($2 + 500) " ""\t"($3 - 150) " ""\t"$4"\t"$5
}' $tss > Parts/tss.SecondPart.bed

awk -F'\t' '{print $1"\t"($2 + 500) " ""\t"($3 - 150) " ""\t"$4"\t"$5
}' $tsP > Parts/tsP.SecondPart.bed

awk -F'\t' '{print $1"\t"($2 + 500) " ""\t"($3 - 150) " ""\t"$4"\t"$5
}' $tsM > Parts/tsM.SecondPart.bed

awk -F'\t' '{print $1"\t"($2 + 350) " ""\t"$3"\t"$4"\t"$5
}' $con > Parts/con.SecondPart.bed

awk -F'\t' '{print $1"\t"($2 + 350) " ""\t"$3"\t"$4"\t"$5
}' $ind > Parts/ind.SecondPart.bed

awk -F'\t' '{print $1"\t"($2 + 350) " ""\t"$3"\t"$4"\t"$5
}' $non > Parts/non.SecondPart.bed



computeMatrix scale-regions -S $pluDht $minDht $pluDms $minDms -R Parts/* -a 0 -b 0 --binSize 350 --transcript_id_designator ACTB --outFileNameMatrix PlusMinusCREMat.tab --outFileSortedRegions PlusMinusCREMatRegions.bed -o mat/groCREPartDensity.mat.gz -p 30 -m 350

plotHeatmap -m mat/groCREPartDensity.mat.gz -out groHeatmapCREParts.pdf --heatmapWidth 10 --heatmapHeight 70 --colorMap "YlGnBu"

cd tmp

tail -n +4 ../PlusMinusCREMat.tab > CREMat.tab
echo -e "groseq_dht.plus\tgroseq_dht.minus\tgroseq_dmso.plus\tgroseq_dmso.minus"  | cat - CREMat.tab > headCREmat.tab
paste ../PlusMinusCREMatRegions.bed headCREmat.tab > concatCRE.tab


chr1    9903362 9903712 chr1:9903362-9903712    .       .       9903362 9903712 0       1       350     9903353 con.FirstPart.bed       3.687   11.81   5.992   6.965
chr1    9903712 9904062 chr1:9903712-9904062    .       .       9903712 9904062 0       1       350     9903703 con.SecondPart.bed      19.84   2.338   27.93   3.847

chr16   2802289 2802639 chr16:2802289 -2802639  .       .       2802289 2802639 0       1       350     2802287 tsP.FirstPart.bed       5.002   13.74   2.492   9.719
chr16   2802639 2802989 chr16:2802639 -2802989  .       .       2802639 2802989 0       1       350     2802637 tsP.SecondPart.bed      1.191   0.1761  0.5544  1.23

chr1    3370879 3371229 chr1:3370879 -3371229   .       .       3370879 3371229 0       1       350     3370876 tsP.FirstPart.bed       0       2.946   0.1806  2.181
chr1    3371229 3371579 chr1:3371229 -3371579   .       .       3371229 3371579 0       1       350     3371226 tsP.SecondPart.bed      2.144   2.135   2.934   0





















# FirstPart
awk -F'\t' '{print $1"\t"($2 + 150) " ""\t"($3 - 500) " ""\t"$4"\t"$5
}' $tsP > Parts/tsP.FirstPart.bed

awk -F'\t' '{print $1"\t"($2 + 150) " ""\t"($3 - 500) " ""\t"$4"\t"$5
}' $tsM > Parts/tsM.FirstPart.bed

awk -F'\t' '{print $1"\t"$2"\t"($3 - 350) " ""\t"$4"\t"$5
}' $coP > Parts/coP.FirstPart.bed

awk -F'\t' '{print $1"\t"$2"\t"($3 - 350) " ""\t"$4"\t"$5
}' $inP > Parts/inP.FirstPart.bed

awk -F'\t' '{print $1"\t"$2"\t"($3 - 350) " ""\t"$4"\t"$5
}' $noP > Parts/noP.FirstPart.bed

awk -F'\t' '{print $1"\t"$2"\t"($3 - 350) " ""\t"$4"\t"$5
}' $coM > Parts/coM.FirstPart.bed

awk -F'\t' '{print $1"\t"$2"\t"($3 - 350) " ""\t"$4"\t"$5
}' $inM > Parts/inM.FirstPart.bed

awk -F'\t' '{print $1"\t"$2"\t"($3 - 350) " ""\t"$4"\t"$5
}' $noM > Parts/noM.FirstPart.bed


# SecondPart
awk -F'\t' '{print $1"\t"($2 + 500) " ""\t"($3 - 150) " ""\t"$4"\t"$5
}' $tsP > Parts/tsP.SecondPart.bed

awk -F'\t' '{print $1"\t"($2 + 500) " ""\t"($3 - 150) " ""\t"$4"\t"$5
}' $tsM > Parts/tsM.SecondPart.bed

awk -F'\t' '{print $1"\t"($2 + 350) " ""\t"$3"\t"$4"\t"$5
}' $coP > Parts/coP.SecondPart.bed

awk -F'\t' '{print $1"\t"($2 + 350) " ""\t"$3"\t"$4"\t"$5
}' $inP > Parts/inP.SecondPart.bed

awk -F'\t' '{print $1"\t"($2 + 350) " ""\t"$3"\t"$4"\t"$5
}' $noP > Parts/noP.SecondPart.bed

awk -F'\t' '{print $1"\t"($2 + 350) " ""\t"$3"\t"$4"\t"$5
}' $coM > Parts/coM.SecondPart.bed

awk -F'\t' '{print $1"\t"($2 + 350) " ""\t"$3"\t"$4"\t"$5
}' $inM > Parts/inM.SecondPart.bed

awk -F'\t' '{print $1"\t"($2 + 350) " ""\t"$3"\t"$4"\t"$5
}' $noM > Parts/noM.SecondPart.bed
