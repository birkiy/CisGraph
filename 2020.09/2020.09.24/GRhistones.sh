




outPwd=/home/ualtintas/GR/histoneSignal



homePwd=/home/ualtintas/GR
cres=(
"GR.con"
"GR.ind"
"GR.non"
"negativeControlGR.final"
"TSS.hg19.Idx.woS"
)

gros=(
"GSM1003453_hg19_wgEncodeBroadHistoneA549H3k04me1Etoh02Sig.bigWig"
"GSM1003495_hg19_wgEncodeBroadHistoneA549H3k04me1Dex100nmSig.bigWig"
"GSM1003578_hg19_wgEncodeBroadHistoneA549H3k27acEtoh02Sig.bigWig"
"GSM1003493_hg19_wgEncodeBroadHistoneA549H3k27acDex100nmSig.bigWig"
"GSM1003561_hg19_wgEncodeBroadHistoneA549H3k04me3Etoh02Sig.bigWig"
"GSM1003542_hg19_wgEncodeBroadHistoneA549H3k04me3Dex100nmSig.bigWig"
"GSM803360_hg19_wgEncodeHaibTfbsA549Pol2Pcr2xEtoh02RawRep1.bigWig"
"GSM803360_hg19_wgEncodeHaibTfbsA549Pol2Pcr2xEtoh02RawRep2.bigWig"
"GSM803361_hg19_wgEncodeHaibTfbsA549Pol2Pcr2xDex100nmRawRep1.bigWig"
"GSM803361_hg19_wgEncodeHaibTfbsA549Pol2Pcr2xDex100nmRawRep2.bigWig"
)

mv SRX186649.bw A549.H3K4me1.Etoh.bw
mv SRX186691.bw A549.H3K4me1.Dex.bw
mv SRX186774.bw A549.H3K27ac.Etoh.bw
mv SRX186689.bw A549.H3K27ac.Dex.bw
mv SRX186757.bw A549.H3K4me3.Etoh.bw
mv SRX186738.bw A549.H3K4me3.Dex.bw
mv SRX100405.bw A549.RNAP2.Etoh.bw
mv SRX100406.bw A549.RNAP2.Dex.bw
mv SRX100418.bw A549.GR.Etoh.bw
mv SRX100416.bw A549.GR.Dex.bw

gros=(
"A549.H3K4me1.Etoh.bw"
"A549.H3K4me1.Dex.bw"
"A549.H3K27ac.Etoh.bw"
"A549.H3K27ac.Dex.bw"
"A549.H3K4me3.Etoh.bw"
"A549.H3K4me3.Dex.bw"
"A549.RNAP2.Etoh.bw"
"A549.RNAP2.Dex.bw"
)

cut -f 1,2,3 $homePwd/negativeControlGR.bed | awk -v s=350 '{print $1, $2-s, $3+s, "negativeControlGR."NR}' | tr ' ' '\t' > $homePwd/negativeControlGR.final.bed


for i in ${cres[@]}
do
  echo $i;
  cut -f4 $homePwd/$i".bed" > "$i.tmp"
  for j in ${gros[@]}
  do
    echo $j;
    bigWigAverageOverBed $homePwd/chipAtlas/$j $homePwd/$i".bed" $outPwd/$i$j".tab";
    echo "bigWigAverageOverBed DONE";
    sort -k 1 $outPwd/$i$j".tab" | cut -f4 | paste -d"\t" "$i.tmp" - > "$i.tab";
    echo "processing DONE";
    cat "$i.tab" > "$i.tmp";
    rm $outPwd/$i$j".tab"
  done
  rm -f "$i.tmp"
done
