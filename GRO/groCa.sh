


# GRO experiments

minDht=/groups/lackgrp/ll_members/common/LNCaP_bigwigs/groseq_analysis/mapping/groseq_dht.minus.bw
pluDht=/groups/lackgrp/ll_members/common/LNCaP_bigwigs/groseq_analysis/mapping/groseq_dht.plus.bw


minDms=/groups/lackgrp/ll_members/common/LNCaP_bigwigs/groseq_analysis/mapping/groseq_dmso.minus.bw
pluDms=/groups/lackgrp/ll_members/common/LNCaP_bigwigs/groseq_analysis/mapping/groseq_dmso.plus.bw


# Promoters
tsP=~/genomeAnnotations/TSS.hg19.+.bed
tsM=~/genomeAnnotations/TSS.hg19.-.bed
pro=~/ARBSs/regions/promoters_ann_5kb.bed

# ARBSs

con=~/ARBSs/regions/cons-arbs.bed
ind=~/ARBSs/regions/ind-arbs.bed
non=~/ARBSs/regions/Non-Active-ARBS.bed


computeMatrix reference-point -S $pluDht $minDht $pluDms $minDms -R $pro $tsP $tsM --referencePoint center -a 3000 -b 3000 -o mat/groPromotersComparison.mat.gz -p 30

plotHeatmap -m mat/groPromotersComparison.mat.gz -out groHeatmapPromoterComparison.pdf --heatmapWidth 10 --heatmapHeight 70 --colorMap "YlGnBu"





computeMatrix reference-point -S $pluDht $minDht $pluDms $minDms -R $tsP $tsM $con $ind $non  --referencePoint center -a 3000 -b 3000 -o mat/groARBSvsPromoters.mat.gz -p 30

plotHeatmap -m mat/groARBSvsPromoters.mat.gz -out groHeatmapARBSvsPromoters.pdf --heatmapWidth 10 --heatmapHeight 70 --colorMap "YlGnBu"



computeMatrix reference-point -S $pluDht $minDht $pluDms $minDms -R $tsP $tsM --referencePoint center -a 1000 -b 1000 -o mat/groPromoters.mat.gz -p 30

plotHeatmap -m mat/groPromoters.mat.gz -out groHeatmapPromoters.pdf --heatmapWidth 10 --heatmapHeight 70 --colorMap "YlGnBu"




computeMatrix reference-point -S $pluDht $minDht $pluDms $minDms -R $con $ind $non  --referencePoint center -a 3000 -b 3000 -o mat/groARBS.mat.gz -p 30

plotHeatmap -m mat/groARBS.mat.gz -out groHeatmapARBS.pdf --heatmapWidth 10 --heatmapHeight 70 --colorMap "YlGnBu"
