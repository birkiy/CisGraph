

GRpath=/home/ualtintas/GR

bigWigMerge \
  $GRpath/STARRbw/GSM3131847_input_dna_1.bw \
  $GRpath/STARRbw/GSM3131848_input_dna_2.bw \
  $GRpath/STARRbw/GSM3131849_input_dna_3.bw \
  $GRpath/STARRbw/GSM3131850_input_dna_4.bw \
  $GRpath/STARRbw/GSM3131851_input_dna_5.bw \
  $GRpath/STARRbw/GSM3131852_input_dna_6.bw \
  $GRpath/STARRbw/GSM3131853_input_dna_7.bw \
  $GRpath/STARRbw/GSM3131854_input_dna_8.bw \
  $GRpath/STARRbw/GSM3131855_input_dna_9.bw \
  $GRpath/STARRbw/GSM3131856_input_dna_10.bw \
  $GRpath/STARRbw/GSM3131857_input_dna_11.bw \
  $GRpath/STARRbw/GSM3131858_input_dna_12.bw \
  $GRpath/STARRbw/input.merged.bedGraph


liftOver $GRpath/STARRbw/input.merged.bedGraph \
  $GRpath/STARRbw/hg38ToHg19.over.chain \
  $GRpath/STARRbw/input.merged.hg19.bedGraph \
  $GRpath/STARRbw/unMapped

bedtools sort -i $GRpath/STARRbw/input.merged.hg19.bedGraph > $GRpath/STARRbw/input.merged.hg19.sorted.bedGraph

bedRemoveOverlap $GRpath/STARRbw/input.merged.hg19.sorted.bedGraph $GRpath/STARRbw/input.merged.hg19.sorted.rmOverlap.bedGraph

bedGraphToBigWig $GRpath/STARRbw/input.merged.hg19.sorted.rmOverlap.bedGraph \
  hg19.chrom.sizes $GRpath/STARRbw/input.merged.hg19.bw


bws=(
"GSM3131822_0hr_1"
"GSM3131823_0hr_2"
"GSM3131824_0hr_3"
"GSM3131825_0hr_4"
"GSM3131826_0hr_5"
"GSM3131827_1hr_1"
"GSM3131828_1hr_2"
"GSM3131829_1hr_3"
"GSM3131830_1hr_4"
"GSM3131831_1hr_5"
"GSM3131832_4hr_1"
"GSM3131833_4hr_2"
"GSM3131834_4hr_3"
"GSM3131835_4hr_4"
"GSM3131836_4hr_5"
"GSM3131837_8hr_1"
"GSM3131838_8hr_2"
"GSM3131839_8hr_3"
"GSM3131840_8hr_4"
"GSM3131841_8hr_5"
"GSM3131842_12hr_1"
"GSM3131843_12hr_2"
"GSM3131844_12hr_3"
"GSM3131845_12hr_4"
"GSM3131846_12hr_5"
)


for bw in ${bws[@]}
do
  echo $bw;
  bigWigToBedGraph $GRpath/STARRbw/$bw".bw" $GRpath/STARRbw/$bw".bedGraph";
  liftOver $GRpath/STARRbw/$bw".bedGraph" \
    $GRpath/STARRbw/hg38ToHg19.over.chain \
    $GRpath/STARRbw/hg19/$bw".hg19.bedGraph" \
    $GRpath/STARRbw/UnMapped/$bw"unMapped";
  bedtools sort -i $GRpath/STARRbw/hg19/$bw".hg19.bedGraph" > $GRpath/STARRbw/hg19/$bw".hg19.sorted.bedGraph";
  bedRemoveOverlap $GRpath/STARRbw/hg19/$bw".hg19.sorted.bedGraph" $GRpath/STARRbw/hg19/$bw".hg19.sorted.rmOverlap.bedGraph";
  bedGraphToBigWig $GRpath/STARRbw/hg19/$bw".hg19.sorted.bedGraph" \
    hg19.chrom.sizes $GRpath/STARRbw/hg19/$bw".hg19.bw";
  echo "lifted";
  bigwigCompare \
    -b1 $GRpath/STARRbw/hg19/$bw".hg19.bw" \
    -b2 $GRpath/STARRbw/input.merged.hg19.bw \
    --outFileName $GRpath/STARRbw/final/$bw".final.bw" \
    --outFileFormat bigwig;
  echo "Done\n"
done




computeMatrix reference-point -S \
  $GRpath/STARRbw/final/"GSM3131822_0hr_1.final.bw" \
  $GRpath/STARRbw/final/"GSM3131823_0hr_2.final.bw" \
  $GRpath/STARRbw/final/"GSM3131824_0hr_3.final.bw" \
  $GRpath/STARRbw/final/"GSM3131825_0hr_4.final.bw" \
  $GRpath/STARRbw/final/"GSM3131826_0hr_5.final.bw" \
  $GRpath/STARRbw/final/"GSM3131827_1hr_1.final.bw" \
  $GRpath/STARRbw/final/"GSM3131828_1hr_2.final.bw" \
  $GRpath/STARRbw/final/"GSM3131829_1hr_3.final.bw" \
  $GRpath/STARRbw/final/"GSM3131830_1hr_4.final.bw" \
  $GRpath/STARRbw/final/"GSM3131831_1hr_5.final.bw" \
  $GRpath/STARRbw/final/"GSM3131832_4hr_1.final.bw" \
  $GRpath/STARRbw/final/"GSM3131833_4hr_2.final.bw" \
  $GRpath/STARRbw/final/"GSM3131834_4hr_3.final.bw" \
  $GRpath/STARRbw/final/"GSM3131835_4hr_4.final.bw" \
  $GRpath/STARRbw/final/"GSM3131836_4hr_5.final.bw" \
  $GRpath/STARRbw/final/"GSM3131837_8hr_1.final.bw" \
  $GRpath/STARRbw/final/"GSM3131838_8hr_2.final.bw" \
  $GRpath/STARRbw/final/"GSM3131839_8hr_3.final.bw" \
  $GRpath/STARRbw/final/"GSM3131840_8hr_4.final.bw" \
  $GRpath/STARRbw/final/"GSM3131841_8hr_5.final.bw" \
  $GRpath/STARRbw/final/"GSM3131842_12hr_1.final.bw" \
  $GRpath/STARRbw/final/"GSM3131843_12hr_2.final.bw" \
  $GRpath/STARRbw/final/"GSM3131844_12hr_3.final.bw" \
  $GRpath/STARRbw/final/"GSM3131845_12hr_4.final.bw" \
  $GRpath/STARRbw/final/"GSM3131846_12hr_5.final.bw" \
  -R $GRpath/GR.con.bed $GRpath/GR.ind.bed $GRpath/GR.non.bed \
  --skipZeros --referencePoint center \
  -a 500 -b 500 -p 30 -q  \
  -o $GRpath/matrices/A549.STARR.mat \
  --outFileNameMatrix $GRpath/matrices/A549.STARR.tab \
  --outFileSortedRegions $GRpath/matrices/A549.STARR.tab



plotHeatmap -m $GRpath/matrices/A549.GRBS.mat \
  -out $GRpath/matrices/A549.GRBS.pdf \
  --colorMap "Purples Purples Purples Purples Purples \
  Blues Blues Blues Blues Blues \
  Greens Greens Greens Greens Greens \
  Oranges Oranges Oranges Oranges Oranges \
  Reds Reds Reds Reds Reds" --missingDataColor 1 \
  --whatToShow "heatmap and colorbar"



















computeMatrix reference-point -S \
  $GRpath/chipAtlas/A549.GR.Dex.bw \
  $GRpath/chipAtlas/A549.GR.Etoh.bw \
  $GRpath/chipAtlas/A549.H3K27ac.Dex.bw \
  $GRpath/chipAtlas/A549.H3K27ac.Etoh.bw \
  $GRpath/chipAtlas/A549.H3K4me1.Dex.bw \
  $GRpath/chipAtlas/A549.H3K4me1.Etoh.bw \
  $GRpath/chipAtlas/A549.H3K4me3.Dex.bw \
  $GRpath/chipAtlas/A549.H3K4me3.Etoh.bw \
  $GRpath/chipAtlas/A549.RNAP2.Dex.bw \
  $GRpath/chipAtlas/A549.RNAP2.Etoh.bw \
  -R $GRpath/GR.con.bed $GRpath/GR.ind.bed $GRpath/GR.non.bed \
  --skipZeros --referencePoint center \
  -a 500 -b 500 -p 30 -q  \
  -o $GRpath/matrices/A549.GRBS.mat \
  --outFileNameMatrix $GRpath/matrices/A549.GRBS.tab \
  --outFileSortedRegions $GRpath/matrices/A549.GRBS.tab


plotHeatmap -m $GRpath/matrices/A549.GRBS.mat -out $GRpath/matrices/A549.GRBS.pdf --colorMap "YlGnBu" --whatToShow "heatmap and colorbar"
