

ssh ualtintas@10.19.148.250

# Working Directory: /home/ualtintas/genomeAnnotations
# Script Location: Cisgraph/Vers2.0/Data/Regions/TSS.sh

# Homer Annotation
annotatePeaks.pl tss hg19 -size -500,500 -cTSS > TSS.hg19.gtf
annotatePeaks.pl tts hg19 -size -500,500 -cTSS > TTS.hg19.gtf

# Convert gtf to BED
cut TSS.hg19.gtf -f2,3,4,16 > tss.pos
# cut TSS.hg19.gtf -f1 > tss.id
cut TSS.hg19.gtf -f5 > tss.st

paste tss.pos tss.st > TSS.hg19.bed
rm tss*
wc -l TSS.hg19.bed

# Remove Haplo Groups
awk '!/hap/{ print $0 }' TSS.hg19.bed > TSS.hg19.woHap.bed
tail -n +2 TSS.hg19.woHap.bed > TSS.hg19.noHead.woHap.bed
bedtools sort -i TSS.hg19.noHead.woHap.bed > TSS.hg19.soted.woHap.bed

# Multiple Indexing
awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4".0\t"$5}' TSS.hg19.soted.woHap.bed > TSS.hg19.soted.woHap.w0.bed


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
} END { for (key in nameDict) {print chrDict[key]"\t"startDict[key]"\t"endDict[key]"\t"key"\t"strDict[key]}}' TSS.hg19.soted.woHap.w0.bed > TSS.hg19.Idx.bed


# Strand Seperation
awk -F'\t' '{if($5 == "+") {print}}' TSS.hg19.Idx.bed > TSS.hg19.+.bed
awk -F'\t' '{if($5 == "-") {print}}' TSS.hg19.Idx.bed > TSS.hg19.-.bed

# Move them to regions folder
mkdir Regions
mv TSS* Regions/.

# Get FASTA
# Working Directory: /home/ualtintas/genomeAnnotations

tss=~/genomeAnnotations/Regions/TSS.hg19.Idx.bed
fi=~/genomeAnnotations/hg19.fa

tsP=~/genomeAnnotations/Regions/TSS.hg19.+.bed
tsM=~/genomeAnnotations/Regions/TSS.hg19.-.bed


bedtools getfasta -fi $fi -bed $tss -fo Fasta/tss.fasta
bedtools getfasta -fi $fi -bed $tsP -fo Fasta/tsP.fasta
bedtools getfasta -fi $fi -bed $tsM -fo Fasta/tsM.fasta
