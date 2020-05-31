

# ssh ualtintas@10.19.148.250

# Working Directory: /home/ualtintas/genomeAnnotations
# Script Location: Cisgraph/Vers2.0/Data/Regions/TTS.sh

# Homer Annotation
annotatePeaks.pl tts hg19 -size -500,500 -cTSS > TTS.hg19.gtf

# Convert gtf to BED
cut TTS.hg19.gtf -f2,3,4,16 > tts.pos
# cut TSS.hg19.gtf -f1 > tss.id
cut TTS.hg19.gtf -f5 > tts.st

paste tts.pos tts.st > TTS.hg19.bed
# rm tts*
wc -l TTS.hg19.bed

# Remove Haplo Groups
awk '!/hap/{ print $0 }' TTS.hg19.bed > TTS.hg19.woHap.bed
tail -n +2 TTS.hg19.woHap.bed > TTS.hg19.noHead.woHap.bed
bedtools sort -i TTS.hg19.noHead.woHap.bed > TTS.hg19.soted.woHap.bed

# Multiple Indexing
awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4".0\t"$5}' TTS.hg19.soted.woHap.bed > TTS.hg19.soted.woHap.w0.bed


awk '{ gene = $4
if (!(gene in nameDict))
{
  nameDict[gene] = 0;
  chrDict[gene] = $1;
  strDict[gene] = $5;
  startDict[gene] = $2;
  endDict[gene] = $3;
}
else
{
  nameDict[gene] += 1
  for(i=length;i!=0;i--)x=x substr(gene,i,1);
  idx=index(x,".");
  idx2=length(gene) - idx + 1;
  for(i=length;i!=0;i--)x=x substr(gene,i,1);
  geneIdx=nameDict[gene];
  geneName=substr(gene,1,idx2-1)"."geneIdx;
  nameDict[geneName] = geneName;
  chrDict[geneName] = $1; strDict[geneName] = $5; startDict[geneName] = $2; endDict[geneName] = $3;
}
} END { for (key in nameDict) {print chrDict[key]"\t"startDict[key]"\t"endDict[key]"\t"key"\t"strDict[key]}}' TTS.hg19.soted.woHap.w0.bed > TTS.hg19.Idx.bed


# Strand Seperation
awk -F'\t' '{if($5 == "+") {print}}' TTS.hg19.Idx.bed > TTS.hg19.+.bed
awk -F'\t' '{if($5 == "-") {print}}' TTS.hg19.Idx.bed > TTS.hg19.-.bed

# Move them to regions folder
# mkdir Regions
mv TTS* Regions/.

# Get FASTA
# Working Directory: /home/ualtintas/genomeAnnotations

tts=~/genomeAnnotations/Regions/TTS.hg19.Idx.bed
fi=~/genomeAnnotations/hg19.fa

ttP=~/genomeAnnotations/Regions/TTS.hg19.+.bed
ttM=~/genomeAnnotations/Regions/TTS.hg19.-.bed


bedtools getfasta -fi $fi -bed $tts -fo Fasta/tts.fasta
bedtools getfasta -fi $fi -bed $ttP -fo Fasta/ttP.fasta
bedtools getfasta -fi $fi -bed $ttM -fo Fasta/ttM.fasta
