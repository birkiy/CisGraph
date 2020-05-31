


chr1	9845	10546


ID=SRX6717234;
Name=NA%20(@%20LNCAP);
Title=ChIP-seq%20of%20RNAP2%20in%20LNCaP/MYC%20KD%20512%20Rep%201;
Cell%20group=Prostate;<br>isolate=LNCaP;
age=0;
biomaterial_provider=Adam%20Sowalsky,%20National%20Cancer%20Institute,%2037%20Convent%20Drive,%20Building%2037%20Room%201062B,%20Bethesda,%20MD%2020892;
sex=male;
tissue=prostate;
cell_line=LNCaP;
treatment=MYC%20knockdown%20construct%201;



1000	.	9845	10546	255,0,0


chr1	9845	10553	ID=SRX6717233;

Name=NA%20(@%20LNCAP);
Title=ChIP-seq%20of%20RNAP2%20in%20LNCaP/Control%20Rep%201;
Cell%20group=Prostate;
<br>isolate=LNCaP;
age=0;
biomaterial_provider=Adam%20Sowalsky,%20National%20Cancer%20Institute,%2037%20Convent%20Drive,%20Building%2037%20Room%201062B,%20Bethesda,%20MD%2020892;
sex=male;
tissue=prostate;
cell_line=LNCaP;
treatment=scrambled%20control;


1000	.	9845	10553	255,0,0


chr1	9999	10247

ID=SRX3474297;
Name=Epitope%20tags%20(@%20LNCAP);
Title=GSM2891162:%20LNCaP%20HA-TRIM24%20%2B%2010nM%20R1881%3B%20Homo%20sapiens%3B%20ChIP-Seq;
Cell%20group=Prostate;
<br>source_name=Prostate%20cancer%20cell%20line;
cell%20line=LNCaP;
genotype=WT;
antibody=HA;
vendor=Abcam;
catalog%20number=ab9110;
lot/batch%20number=GR261166-6;

1000	.	9999	10247	255,0,0


###########################################
#                                         #
#          Umut Berkay Altıntaş           #
#                                         #
###########################################


outPath=/kuacc/users/ualtintas20/LNCaP/BED/
lnBed=/kuacc/users/ualtintas20/LNCaP/ALL.Prs.50.AllAg.LNCAP.bed
lnBed=/kuacc/users/ualtintas20/LNCaP/ALL.Prs.05.AllAg.LNCAP.bed
# Histones=$@


Histones=("H3K27Ac" "H3K4me3" "H3K4me2" "H3K4me1" "H3ac" "H4ac" "H3K36me3" "H2AZ" "H2AZac" "H3" "H3K9me3")

for histone in ${Histones[@]} ; do
  histoneS=$histone"%20";
  echo $histoneS $histone
	grep $histoneS $lnBed | awk -F'\t' -v var="$histone" '{print $1"\t"$2"\t"$3"\t"var"."NR }' > $outPath"Histones/"$histone".bed"
done

grep "H3K27Ac%20" $lnBed
grep "H3K4me3%20" $lnBed
grep "H3K4me2%20" $lnBed
grep "H3K4me1%20" $lnBed
grep "H3ac%20" $lnBed
grep "H4ac%20" $lnBed
grep "H3K36me3%20" $lnBed
grep "H2AZ%20" $lnBed
grep "H2AZac%20" $lnBed
grep "H3%20" $lnBed
grep "H3K9me3%20" $lnBed


Proteins=("RNAP2" "AR" "FOXA1" "ARID1A" "CHD1" "CHD4" "CSNK2A1" "CTBP1" "CTBP2" "EZH2" "GATA2" "GRHL2" "HOXB13" "KDM1A" "MTOR" "NMYC" "PIAS1" "POU2F1" "RELA" "SFPQ" "SMARCA1" "SMARCA2" "SMARCA4" "SMARCA5" "SUZ12" "TCF7L2" "TET2" "CTCF" "TLE3" "TRIM24" "TRIM28" "WDHD1" "WDR5")

for protein in ${Proteins[@]} ; do
  proteinS=$protein"%20";
  echo $proteinS $protein
	grep $proteinS $lnBed | awk -F'\t' -v var="$protein" '{print $1"\t"$2"\t"$3"\t"var"."NR }' > $outPath"Proteins/"$protein".bed";
done



protein=HIF1A
proteinS="HIF1A%20"
echo $proteinS $protein
grep $proteinS $lnBed | awk -F'\t' -v var="$protein" '{print $1"\t"$2"\t"$3"\t"var"."NR }' > $outPath"Proteins/"$protein".bed";

protein=MED
proteinS="MED"
echo $proteinS $protein
grep $proteinS $lnBed | awk -F'\t' -v var="$protein" '{print $1"\t"$2"\t"$3"\t"var"."NR }' > $outPath"Proteins/"$protein".bed";

protein=NKX3
proteinS="NKX3.1%20"
echo $proteinS $protein
grep $proteinS $lnBed | awk -F'\t' -v var="$protein" '{print $1"\t"$2"\t"$3"\t"var"."NR }' > $outPath"Proteins/"$protein".bed";

protein=TFAP4
proteinS="TFAP4%20"
echo $proteinS $protein
grep $proteinS $lnBed | awk -F'\t' -v var="$protein" '{print $1"\t"$2"\t"$3"\t"var"."NR }' > $outPath"Proteins/"$protein".bed";






grep "RNAP2%20" $lnBed
grep "AR%20" $lnBed
grep "FOXA1%20" $lnBed
grep "ARID1A%20" $lnBed
grep "CHD1%20" $lnBed
grep "CHD4%20" $lnBed
grep "CSNK2A1%20" $lnBed
grep "CTBP1%20" $lnBed
grep "CTBP2%20" $lnBed
grep "EZH2%20" $lnBed
grep "GATA2%20" $lnBed
grep "GRHL2%20" $lnBed
grep "HOXB13%20" $lnBed
grep "KDM1A%20" $lnBed
grep "MRE11A%20" $lnBed
grep "MTOR%20" $lnBed
grep "NMYC%20" $lnBed
grep "PIAS1%20" $lnBed
grep "POU2F1%20" $lnBed
grep "RELA%20" $lnBed
grep "SFPQ%20" $lnBed
grep "SMARCA1%20" $lnBed
grep "SMARCA2%20" $lnBed
grep "SMARCA4%20" $lnBed
grep "SMARCA5%20" $lnBed
grep "SUZ12%20" $lnBed
grep "TCF7L2%20" $lnBed
grep "TET2%20" $lnBed
grep "CTCF%20" $lnBed
grep "TLE3%20" $lnBed
grep "TRIM24%20" $lnBed
grep "TRIM28%20" $lnBed
grep "WDHD1%20" $lnBed
grep "WDR5%20" $lnBed



grep "HIF1A%20" $lnBed > BED/HIF1A.bed
grep "MED" $lnBed > BED/MED.bed
grep "NKX3.1%20" $lnBed > BED/NKX31.bed
grep "TFAP4%20" $lnBed > BED/TFAP4.bed


awk -F'\t' ''








#  These are not exist
# grep "FOXP%20" $lnBed
# grep "GSE62492" $lnBed
# grep "RUNX1%20" $lnBed
# gerp "GSE62492" $lnBed
#
# grep "MRE11A%20" $lnBed
# grep "GSE63202" $lnBed
