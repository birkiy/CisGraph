

minDht=/groups/lackgrp/ll_members/common/LNCaP_bigwigs/groseq_analysis/mapping/groseq_dht.minus.bw
pluDht=/groups/lackgrp/ll_members/common/LNCaP_bigwigs/groseq_analysis/mapping/groseq_dht.plus.bw


minDms=/groups/lackgrp/ll_members/common/LNCaP_bigwigs/groseq_analysis/mapping/groseq_dmso.minus.bw
pluDms=/groups/lackgrp/ll_members/common/LNCaP_bigwigs/groseq_analysis/mapping/groseq_dmso.plus.bw

# Promoters
tsP=~/genomeAnnotations/Regions/TSS.hg19.+.bed
tsM=~/genomeAnnotations/Regions/TSS.hg19.-.bed

# FirstPart
awk -F'\t' '{print $1"\t"$2"\t"($3 - 500) " ""\t"$4"\t"$5
}' $tsP > tsP.FirstPart.bed

tsPF=~/genomeAnnotations/Regions/tsP.FirstPart.bed

awk -F'\t' '{print $1"\t"$2"\t"($3 - 500) " ""\t"$4"\t"$5
}' $tsM > tsM.FirstPart.bed

tsMF=~/genomeAnnotations/Regions/tsM.FirstPart.bed

# SecondPart
awk -F'\t' '{print $1"\t"($2 + 500) " ""\t"$3"\t"$4"\t"$5
}' $tsP > tsP.SecondPart.bed

tsPS=~/genomeAnnotations/Regions/tsP.SecondPart.bed
OvohF5el

awk -F'\t' '{print $1"\t"($2 + 500) " ""\t"$3"\t"$4"\t"$5
}' $tsM > tsM.SecondPart.bed

tsMS=~/genomeAnnotations/Regions/tsM.SecondPart.bed


# #tsP.FirstPart.bed:27347	tsM.FirstPart.bed:26083	tsP.SecondPart.bed:27375	tsM.SecondPart.bed:25993
# #downstream:0	upstream:0	body:1000	bin size:500	unscaled 5 prime:0	unscaled 3 prime:0
# tsP.FirstPart.bed:27347	tsM.FirstPart.bed:26083	tsP.SecondPart.bed:27375	tsM.SecondPart.bed:25993



computeMatrix scale-regions -S $pluDht $minDht $pluDms $minDms -R $tsPF $tsMF $tsPS $tsMS -a 0 -b 0 --binSize 500 --transcript_id_designator ACTB --outFileNameMatrix PlusMinusMat.tab --outFileSortedRegions PlusMinusMatRegions.bed -o mat/groPromotersPartDensity.mat.gz -p 30 --sortRegions descend --skipZeros


plotHeatmap -m mat/groPromotersPartDensity.mat.gz -out groHeatmapPromoterParts.pdf --heatmapWidth 10 --heatmapHeight 70 --colorMap "YlGnBu"



conMF=~/ARBSs/regions/con.ARBS.M.bed
conPS=~/ARBSs/regions/con.ARBS.P.bed

indMF=~/ARBSs/regions/ind.ARBS.M.bed
indPS=~/ARBSs/regions/ind.ARBS.P.bed

nonMF=~/ARBSs/regions/non.ARBS.M.bed
nonPS=~/ARBSs/regions/non.ARBS.P.bed

computeMatrix scale-regions -S $pluDht $minDht $pluDms $minDms -R $conMF $conPS $indMF $indPS $nonMF $nonPS -a 0 -b 0 --binSize 350 --transcript_id_designator ACTB --outFileNameMatrix PlusMinusArbsMat.tab --outFileSortedRegions PlusMinusArbsMatRegions.bed -o mat/groARBSPartDensity.mat.gz -p 30 --sortRegions descend --skipZeros -m 700

# computeMatrix scale-regions -S $pluDht $minDht $pluDms $minDms -R $arbPS $arbMF -a 0 -b 0 --binSize 350 --transcript_id_designator ACTB --outFileNameMatrix PlusMinusArbsMat.tab --outFileSortedRegions PlusMinusArbsMatRegions.bed -o mat/groARBSPartDensity.mat.gz -p 30 --sortRegions descend --skipZeros


plotHeatmap -m mat/groARBSPartDensity.mat.gz -out groHeatmapPromoterParts.pdf --heatmapWidth 10 --heatmapHeight 70 --colorMap "YlGnBu"
