

TFpath=/kuacc/users/ualtintas20/LNCaP/BED/Proteins


TFs=(AR.bed ARID1A.bed CHD1.bed CHD4.bed CSNK2A1.bed CTBP1.bed CTBP2.bed CTCF.bed EZH2.bed FOXA1.bed GATA2.bed GRHL2.bed HIF1A.bed HOXB13.bed KDM1A.bed MED.bed MTOR.bed NKX31.bed NKX3.bed NMYC.bed PIAS1.bed POU2F1.bed RELA.bed RNAP2.bed SFPQ.bed SMARCA1.bed SMARCA2.bed SMARCA4.bed SMARCA5.bed SUZ12.bed TCF7L2.bed TET2.bed TFAP4.bed TLE3.bed TRIM24.bed TRIM28.bed WDHD1.bed WDR5.bed)

for tf in ${TFs[@]}
do
  bedtools sort -i $TFpath/$tf > $TFpath/processed/sorted.$tf
  bedtools merge -i $TFpath/processed/sorted.$tf > $TFpath/processed/merged.$tf
done


multiIntersectBed -i merged.AR.bed \
merged.ARID1A.bed \
merged.CHD1.bed \
merged.CHD4.bed \
merged.CSNK2A1.bed \
merged.CTBP1.bed \
merged.CTBP2.bed \
merged.CTCF.bed \
merged.EZH2.bed \
merged.FOXA1.bed \
merged.GATA2.bed \
merged.GRHL2.bed \
merged.HIF1A.bed \
merged.HOXB13.bed \
merged.KDM1A.bed \
merged.MED.bed \
merged.MTOR.bed \
merged.NKX31.bed \
merged.NKX3.bed \
merged.NMYC.bed \
merged.PIAS1.bed \
merged.POU2F1.bed \
merged.RELA.bed \
merged.RNAP2.bed \
merged.SFPQ.bed \
merged.SMARCA1.bed \
merged.SMARCA2.bed \
merged.SMARCA4.bed \
merged.SMARCA5.bed \
merged.SUZ12.bed \
merged.TCF7L2.bed \
merged.TET2.bed \
merged.TFAP4.bed \
merged.TLE3.bed \
merged.TRIM24.bed \
merged.TRIM28.bed \
merged.WDHD1.bed \
merged.WDR5.bed > LNCaP.CRE.bed


bedtools intersect -a LNCaP.CRE.bed -b ~/ARBSs/regions/cons-arbs.bed -wa -wb > cre/con.cre.bed
bedtools intersect -a LNCaP.CRE.bed -b ~/ARBSs/regions/ind-arbs.bed -wa -wb > cre/ind.cre.bed
bedtools intersect -a LNCaP.CRE.bed -b ~/ARBSs/regions/Non-Active-ARBS.bed -wa -wb > cre/non.cre.bed
bedtools intersect -a LNCaP.CRE.bed -b ~/ARBSs/regions/negativeControl.ARBS.bed -wa -wb > cre/nAR.cre.bed
bedtools intersect -a LNCaP.CRE.bed -b ~/genomeAnnotations/Regions/TSS.hg19.Idx.bed -wa -wb > cre/tss.cre.bed
