

cut -f1,2,3 $idrPath/22RV1.AR-C19.idr > $idrPath/22RV1.AR-C19.idr.peaks.bed
awk -F'\t' '{print $0"\t22RV1.AR-C19.idr."NR}' $idrPath/22RV1.AR-C19.idr.peaks.bed > $idrPath/22RV1.AR-C19.idr.peaks.bed.tmp
mv $idrPath/22RV1.AR-C19.idr.peaks.bed.tmp $idrPath/22RV1.AR-C19.idr.peaks.bed

cut -f1,2,3 $idrPath/22RV1.AR-V7.idr > $idrPath/22RV1.AR-V7.idr.peaks.bed
awk -F'\t' '{print $0"\t22RV1.AR-V7.idr."NR}' $idrPath/22RV1.AR-V7.idr.peaks.bed > $idrPath/22RV1.AR-V7.idr.peaks.bed.tmp
mv $idrPath/22RV1.AR-V7.idr.peaks.bed.tmp $idrPath/22RV1.AR-V7.idr.peaks.bed

cut -f1,2,3 $idrPath/LN95.AR-C19.idr > $idrPath/LN95.AR-C19.idr.peaks.bed
awk -F'\t' '{print $0"\tLN95.AR-C19.idr."NR}' $idrPath/LN95.AR-C19.idr.peaks.bed > $idrPath/LN95.AR-C19.idr.peaks.bed.tmp
mv $idrPath/LN95.AR-C19.idr.peaks.bed.tmp $idrPath/LN95.AR-C19.idr.peaks.bed

cut -f1,2,3 $idrPath/LN95.AR-V7.idr > $idrPath/LN95.AR-V7.idr.peaks.bed
awk -F'\t' '{print $0"\tLN95.AR-V7.idr."NR}' $idrPath/LN95.AR-V7.idr.peaks.bed > $idrPath/LN95.AR-V7.idr.peaks.bed.tmp
mv $idrPath/LN95.AR-V7.idr.peaks.bed.tmp $idrPath/LN95.AR-V7.idr.peaks.bed

cut -f1,2,3 $idrPath/LNCaP.dht.AR.idr > $idrPath/LNCaP.dht.AR.idr.peaks.bed
awk -F'\t' '{print $0"\tLNCaP.dht.AR.idr."NR}' $idrPath/LNCaP.dht.AR.idr.peaks.bed > $idrPath/LNCaP.dht.AR.idr.peaks.bed.tmp
mv $idrPath/LNCaP.dht.AR.idr.peaks.bed.tmp $idrPath/LNCaP.dht.AR.idr.peaks.bed

cut -f1,2,3 $idrPath/LNCaP.veh.AR.idr > $idrPath/LNCaP.veh.AR.idr.peaks.bed
awk -F'\t' '{print $0"\tLNCaP.veh.AR.idr."NR}' $idrPath/LNCaP.veh.AR.idr.peaks.bed > $idrPath/LNCaP.veh.AR.idr.peaks.bed.tmp
mv $idrPath/LNCaP.veh.AR.idr.peaks.bed.tmp $idrPath/LNCaP.veh.AR.idr.peaks.bed

cut -f1,2,3 $idrPath/malignant.1.AR.idr > $idrPath/malignant.1.AR.idr.peaks.bed
awk -F'\t' '{print $0"\tmalignant.1.AR.idr."NR}' $idrPath/malignant.1.AR.idr.peaks.bed > $idrPath/malignant.1.AR.idr.peaks.bed.tmp
mv $idrPath/malignant.1.AR.idr.peaks.bed.tmp $idrPath/malignant.1.AR.idr.peaks.bed

cut -f1,2,3 $idrPath/malignant.2.AR.idr > $idrPath/malignant.2.AR.idr.peaks.bed
awk -F'\t' '{print $0"\tmalignant.2.AR.idr."NR}' $idrPath/malignant.2.AR.idr.peaks.bed > $idrPath/malignant.2.AR.idr.peaks.bed.tmp
mv $idrPath/malignant.2.AR.idr.peaks.bed.tmp $idrPath/malignant.2.AR.idr.peaks.bed

cut -f1,2,3 $idrPath/malignant.3.AR.idr > $idrPath/malignant.3.AR.idr.peaks.bed
awk -F'\t' '{print $0"\tmalignant.3.AR.idr."NR}' $idrPath/malignant.3.AR.idr.peaks.bed > $idrPath/malignant.3.AR.idr.peaks.bed.tmp
mv $idrPath/malignant.3.AR.idr.peaks.bed.tmp $idrPath/malignant.3.AR.idr.peaks.bed

cut -f1,2,3 $idrPath/malignant.4.AR.idr > $idrPath/malignant.4.AR.idr.peaks.bed
awk -F'\t' '{print $0"\tmalignant.4.AR.idr."NR}' $idrPath/malignant.4.AR.idr.peaks.bed > $idrPath/malignant.4.AR.idr.peaks.bed.tmp
mv $idrPath/malignant.4.AR.idr.peaks.bed.tmp $idrPath/malignant.4.AR.idr.peaks.bed

cut -f1,2,3 $idrPath/non-malignant.1.AR.idr > $idrPath/non-malignant.1.AR.idr.peaks.bed
awk -F'\t' '{print $0"\tnon-malignant.1.AR.idr."NR}' $idrPath/non-malignant.1.AR.idr.peaks.bed > $idrPath/non-malignant.1.AR.idr.peaks.bed.tmp
mv $idrPath/non-malignant.1.AR.idr.peaks.bed.tmp $idrPath/non-malignant.1.AR.idr.peaks.bed

cut -f1,2,3 $idrPath/non-malignant.2.AR.idr > $idrPath/non-malignant.2.AR.idr.peaks.bed
awk -F'\t' '{print $0"\tnon-malignant.2.AR.idr."NR}' $idrPath/non-malignant.2.AR.idr.peaks.bed > $idrPath/non-malignant.2.AR.idr.peaks.bed.tmp
mv $idrPath/non-malignant.2.AR.idr.peaks.bed.tmp $idrPath/non-malignant.2.AR.idr.peaks.bed
