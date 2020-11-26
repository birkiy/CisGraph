samtools merge -@ 20 -r $home/mapping/processed/GR.0h.merged.final.bam \
  $home/mapping/processed/GR.0h.dex.rep1.final.bam \
  $home/mapping/processed/GR.0h.dex.rep2.final.bam \
  $home/mapping/processed/GR.0h.dex.rep3.final.bam \
  $home/mapping/processed/GR.0h.dex.rep4.final.bam \
  $home/mapping/processed/GR.0h.dex.rep5.final.bam


samtools merge -@ 20 -r $home/mapping/processed/GR.1h.merged.final.bam \
  $home/mapping/processed/GR.1h.dex.rep1.final.bam \
  $home/mapping/processed/GR.1h.dex.rep2.final.bam \
  $home/mapping/processed/GR.1h.dex.rep3.final.bam \
  $home/mapping/processed/GR.1h.dex.rep4.final.bam \
  $home/mapping/processed/GR.1h.dex.rep5.final.bam


samtools merge -@ 20 -r $home/mapping/processed/GR.4h.merged.final.bam \
  $home/mapping/processed/GR.4h.dex.rep1.final.bam \
  $home/mapping/processed/GR.4h.dex.rep2.final.bam \
  $home/mapping/processed/GR.4h.dex.rep3.final.bam \
  $home/mapping/processed/GR.4h.dex.rep4.final.bam \
  $home/mapping/processed/GR.4h.dex.rep5.final.bam

samtools merge -@ 20 -r $home/mapping/processed/GR.8h.merged.final.bam \
  $home/mapping/processed/GR.8h.dex.rep1.final.bam \
  $home/mapping/processed/GR.8h.dex.rep2.final.bam \
  $home/mapping/processed/GR.8h.dex.rep3.final.bam \
  $home/mapping/processed/GR.8h.dex.rep4.final.bam \
  $home/mapping/processed/GR.8h.dex.rep5.final.bam


samtools merge -@ 20 -r $home/mapping/processed/GR.12h.merged.final.bam \
  $home/mapping/processed/GR.12h.dex.rep1.final.bam \
  $home/mapping/processed/GR.12h.dex.rep2.final.bam \
  $home/mapping/processed/GR.12h.dex.rep3.final.bam \
  $home/mapping/processed/GR.12h.dex.rep4.final.bam \
  $home/mapping/processed/GR.12h.dex.rep5.final.bam

samtools merge -@ 20 -r $home/mapping/processed/input.GSE114063.merged.final.bam \
  $home/mapping/processed/input.GSE114063.pool10.final.bam \
  $home/mapping/processed/input.GSE114063.pool11.final.bam \
  $home/mapping/processed/input.GSE114063.pool12.final.bam \
  $home/mapping/processed/input.GSE114063.pool1.final.bam \
  $home/mapping/processed/input.GSE114063.pool2.final.bam \
  $home/mapping/processed/input.GSE114063.pool3.final.bam \
  $home/mapping/processed/input.GSE114063.pool4.final.bam \
  $home/mapping/processed/input.GSE114063.pool5.final.bam \
  $home/mapping/processed/input.GSE114063.pool6.final.bam \
  $home/mapping/processed/input.GSE114063.pool7.final.bam \
  $home/mapping/processed/input.GSE114063.pool8.final.bam \
  $home/mapping/processed/input.GSE114063.pool9.final.bam




samtools index $home/mapping/processed/GR.0h.dex.merged.final.bam $home/mapping/processed/GR.0h.dex.merged.final.bam.bai
samtools index $home/mapping/processed/GR.1h.dex.merged.final.bam $home/mapping/processed/GR.1h.dex.merged.final.bam.bai
samtools index $home/mapping/processed/GR.4h.dex.merged.final.bam $home/mapping/processed/GR.4h.dex.merged.final.bam.bai
samtools index $home/mapping/processed/GR.8h.dex.merged.final.bam $home/mapping/processed/GR.8h.dex.merged.final.bam.bai
samtools index $home/mapping/processed/GR.12h.dex.merged.final.bam $home/mapping/processed/GR.12h.dex.merged.final.bam.bai
samtools index $home/mapping/processed/input.GSE114063.merged.final.bam $home/mapping/processed/input.GSE114063.merged.final.bam.bai
