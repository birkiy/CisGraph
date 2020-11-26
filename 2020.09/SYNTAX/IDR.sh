

home=/groups/lackgrp/ll_members/berkay/SYNTAX/analysis/peaks
sort -k8,8nr $home/AR.dht.rep1_peaks.narrowPeak > $home/idr/AR.dht.rep1.sorted.peaks.narrowPeak
sort -k8,8nr $home/AR.dht.rep2_peaks.narrowPeak > $home/idr/AR.dht.rep2.sorted.peaks.narrowPeak
sort -k8,8nr $home/AR.dmso.rep1_peaks.narrowPeak > $home/idr/AR.dmso.rep1.sorted.peaks.narrowPeak
sort -k8,8nr $home/AR.dmso.rep2_peaks.narrowPeak > $home/idr/AR.dmso.rep2.sorted.peaks.narrowPeak
sort -k8,8nr $home/FOXA1.dht.rep1_peaks.narrowPeak > $home/idr/FOXA1.dht.rep1.sorted.peaks.narrowPeak
sort -k8,8nr $home/FOXA1.dht.rep2_peaks.narrowPeak > $home/idr/FOXA1.dht.rep2.sorted.peaks.narrowPeak
sort -k8,8nr $home/FOXA1.dmso.rep1_peaks.narrowPeak > $home/idr/FOXA1.dmso.rep1.sorted.peaks.narrowPeak
sort -k8,8nr $home/FOXA1.dmso.rep2_peaks.narrowPeak > $home/idr/FOXA1.dmso.rep2.sorted.peaks.narrowPeak






idr --samples $home/idr/FOXA1.dht.rep1.sorted.peaks.narrowPeak $home/idr/FOXA1.dht.rep2.sorted.peaks.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file $home/idr/FOXA1.dht.idr \
--plot \
--log-output-file $home/idr/FOXA1.dht.idr.log




idr --samples $home/idr/FOXA1.dmso.rep1.sorted.peaks.narrowPeak $home/idr/FOXA1.dmso.rep2.sorted.peaks.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file $home/idr/FOXA1.dmso.idr \
--plot \
--log-output-file $home/idr/FOXA1.dmso.idr.log







idr --samples $home/idr/AR.dht.rep1.sorted.peaks.narrowPeak $home/idr/AR.dht.rep2.sorted.peaks.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file $home/idr/AR.dht.idr \
--plot \
--log-output-file $home/idr/AR.dht.idr.log



idr --samples $home/idr/AR.dmso.rep1.sorted.peaks.narrowPeak $home/idr/AR.dmso.rep2.sorted.peaks.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file $home/idr/AR.dmso.idr \
--plot \
--log-output-file $home/idr/AR.dmso.idr.log
