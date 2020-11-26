

bigWigToBedGraph links/A549.H3K27ac.dex.8h.rep2.bigWig converted/A549.H3K27ac.dex.8h.rep2.38.bedGraph
liftOver converted/A549.H3K27ac.dex.8h.rep2.38.bedGraph                     /groups/lackgrp/ll_members/berkay/STARRbegin/A549/code/hg38ToHg19.over.chain                     converted/A549.H3K27ac.dex.8h.rep2.19.bedGraph                     converted/A549.H3K27ac.dex.8h.rep2.un.bedGraph
bedtools sort -i converted/A549.H3K27ac.dex.8h.rep2.19.bedGraph > converted/A549.H3K27ac.dex.8h.rep2.so.bedGraph
bedRemoveOverlap converted/A549.H3K27ac.dex.8h.rep2.so.bedGraph converted/A549.H3K27ac.dex.8h.rep2.rm.bedGraph
bedGraphToBigWig converted/A549.H3K27ac.dex.8h.rep2.rm.bedGraph                     /home/ualtintas/genomeAnnotations/hg19.chrom.sizes converted/A549.H3K27ac.dex.8h.rep2.bigWig


bigWigToBedGraph links/A549.H3K27ac.dex.0h.rep1.bigWig converted/A549.H3K27ac.dex.0h.rep1.38.bedGraph
liftOver converted/A549.H3K27ac.dex.0h.rep1.38.bedGraph                     /groups/lackgrp/ll_members/berkay/STARRbegin/A549/code/hg38ToHg19.over.chain                     converted/A549.H3K27ac.dex.0h.rep1.19.bedGraph                     converted/A549.H3K27ac.dex.0h.rep1.un.bedGraph
bedtools sort -i converted/A549.H3K27ac.dex.0h.rep1.19.bedGraph > converted/A549.H3K27ac.dex.0h.rep1.so.bedGraph
bedRemoveOverlap converted/A549.H3K27ac.dex.0h.rep1.so.bedGraph converted/A549.H3K27ac.dex.0h.rep1.rm.bedGraph
bedGraphToBigWig converted/A549.H3K27ac.dex.0h.rep1.rm.bedGraph                     /home/ualtintas/genomeAnnotations/hg19.chrom.sizes converted/A549.H3K27ac.dex.0h.rep1.bigWig



bigWigToBedGraph links/A549.H3K27ac.dex.0h.rep3.bigWig converted/A549.H3K27ac.dex.0h.rep3.38.bedGraph
liftOver converted/A549.H3K27ac.dex.0h.rep3.38.bedGraph                     /groups/lackgrp/ll_members/berkay/STARRbegin/A549/code/hg38ToHg19.over.chain                     converted/A549.H3K27ac.dex.0h.rep3.19.bedGraph                     converted/A549.H3K27ac.dex.0h.rep3.un.bedGraph
bedtools sort -i converted/A549.H3K27ac.dex.0h.rep3.19.bedGraph > converted/A549.H3K27ac.dex.0h.rep3.so.bedGraph
bedRemoveOverlap converted/A549.H3K27ac.dex.0h.rep3.so.bedGraph converted/A549.H3K27ac.dex.0h.rep3.rm.bedGraph
bedGraphToBigWig converted/A549.H3K27ac.dex.0h.rep3.rm.bedGraph                     /home/ualtintas/genomeAnnotations/hg19.chrom.sizes converted/A549.H3K27ac.dex.0h.rep3.bigWig



bigWigAverageOverBed converted/A549.H3K27ac.dex.0h.rep1.bigWig <(cut -f1,2,3,4 /groups/lackgrp/ll_members/berkay/STARRbegin/peaks/con.0.5.GR.8h.bed) coverage/A549.H3K27ac.dex.0h.rep1.con.avg.tab
