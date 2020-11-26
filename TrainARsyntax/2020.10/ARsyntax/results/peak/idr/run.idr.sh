#!/bin/bash
#SBATCH --job-name=bpnet.idr
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=2
#SBATCH --export=all
#SBATCH -p long


idrPath=/groups/lackgrp/ll_members/berkay/ARsyntax/results/peak/idr

idr --samples \
  $idrPath/22RV1.AR-C19.rep1.sorted.peaks.narrowPeak \
  $idrPath/22RV1.AR-C19.rep2.sorted.peaks.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value \
  --output-file $idrPath/22RV1.AR-C19.idr \
  --plot



idr --samples \
  $idrPath/22RV1.AR-V7.rep1.sorted.peaks.narrowPeak \
  $idrPath/22RV1.AR-V7.rep2.sorted.peaks.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value \
  --output-file $idrPath/22RV1.AR-V7.idr \
  --plot



idr --samples \
  $idrPath/LN95.AR-C19.rep1.sorted.peaks.narrowPeak \
  $idrPath/LN95.AR-C19.rep2.sorted.peaks.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value \
  --output-file $idrPath/LN95.AR-C19.idr \
  --plot



idr --samples \
  $idrPath/LN95.AR-V7.rep1.sorted.peaks.narrowPeak \
  $idrPath/LN95.AR-V7.rep2.sorted.peaks.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value \
  --output-file $idrPath/LN95.AR-V7.idr \
  --plot



idr --samples \
  $idrPath/LNCaP.dht.AR.rep1.sorted.peaks.narrowPeak \
  $idrPath/LNCaP.dht.AR.rep2.sorted.peaks.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value \
  --output-file $idrPath/LNCaP.dht.AR.idr \
  --plot


idr --samples \
  $idrPath/LNCaP.veh.AR.rep1.sorted.peaks.narrowPeak \
  $idrPath/LNCaP.veh.AR.rep2.sorted.peaks.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value \
  --output-file $idrPath/LNCaP.veh.AR.idr \
  --plot


idr --samples \
  $idrPath/malignant.1.AR.rep1.sorted.peaks.narrowPeak \
  $idrPath/malignant.1.AR.rep2.sorted.peaks.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value \
  --output-file $idrPath/malignant.1.AR.idr \
  --plot


idr --samples \
  $idrPath/malignant.2.AR.rep1.sorted.peaks.narrowPeak \
  $idrPath/malignant.2.AR.rep2.sorted.peaks.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value \
  --output-file $idrPath/malignant.2.AR.idr \
  --plot

idr --samples \
  $idrPath/malignant.3.AR.rep1.sorted.peaks.narrowPeak \
  $idrPath/malignant.3.AR.rep2.sorted.peaks.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value \
  --output-file $idrPath/malignant.3.AR.idr \
  --plot


idr --samples \
  $idrPath/malignant.4.AR.rep1.sorted.peaks.narrowPeak \
  $idrPath/malignant.4.AR.rep2.sorted.peaks.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value \
  --output-file $idrPath/malignant.4.AR.idr \
  --plot



idr --samples \
  $idrPath/non-malignant.1.AR.rep1.sorted.peaks.narrowPeak \
  $idrPath/non-malignant.1.AR.rep2.sorted.peaks.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value \
  --output-file $idrPath/non-malignant.1.AR.idr \
  --plot


idr --samples \
  $idrPath/non-malignant.2.AR.rep1.sorted.peaks.narrowPeak \
  $idrPath/non-malignant.2.AR.rep2.sorted.peaks.narrowPeak \
  --input-file-type narrowPeak \
  --rank p.value \
  --output-file $idrPath/non-malignant.2.AR.idr \
  --plot
