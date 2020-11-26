#!/bin/bash
#SBATCH --job-name=bpnet.SumBW
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=40
#SBATCH --export=all
#SBATCH -p long



bigwigCompare -b1 22RV1.AR-C19.rep1.+.5end.bigWig \
-b2 22RV1.AR-C19.rep2.+.5end.bigWig \
--operation add --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed -p 40 \
-bs=1 -o 22RV1.AR-C19.+.5end.bigWig \

bigwigCompare -b1 22RV1.AR-C19.rep1.-.5end.bigWig \
-b2 22RV1.AR-C19.rep2.-.5end.bigWig \
--operation add --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed -p 40 \
-bs=1 -o 22RV1.AR-C19.-.5end.bigWig \

bigwigCompare -b1 22RV1.AR-V7.rep1.+.5end.bigWig \
-b2 22RV1.AR-V7.rep2.+.5end.bigWig \
--operation add --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed -p 40 \
-bs=1 -o 22RV1.AR-V7.+.5end.bigWig \

bigwigCompare -b1 22RV1.AR-V7.rep1.-.5end.bigWig \
-b2 22RV1.AR-V7.rep2.-.5end.bigWig \
--operation add --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed -p 40 \
-bs=1 -o 22RV1.AR-V7.-.5end.bigWig \

bigwigCompare -b1 LN95.AR-C19.rep1.+.5end.bigWig \
-b2 LN95.AR-C19.rep2.+.5end.bigWig \
--operation add --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed -p 40 \
-bs=1 -o LN95.AR-C19.+.5end.bigWig \

bigwigCompare -b1 LN95.AR-C19.rep1.-.5end.bigWig \
-b2 LN95.AR-C19.rep2.-.5end.bigWig \
--operation add --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed -p 40 \
-bs=1 -o LN95.AR-C19.-.5end.bigWig \

bigwigCompare -b1 LN95.AR-V7.rep1.+.5end.bigWig \
-b2 LN95.AR-V7.rep2.+.5end.bigWig \
--operation add --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed -p 40 \
-bs=1 -o LN95.AR-V7.+.5end.bigWig \

bigwigCompare -b1 LN95.AR-V7.rep1.-.5end.bigWig \
-b2 LN95.AR-V7.rep2.-.5end.bigWig \
--operation add --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed -p 40 \
-bs=1 -o LN95.AR-V7.-.5end.bigWig \

bigwigCompare -b1 LNCaP.dht.AR.rep1.+.5end.bigWig \
-b2 LNCaP.dht.AR.rep2.+.5end.bigWig \
--operation add --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed -p 40 \
-bs=1 -o LNCaP.dht.AR.+.5end.bigWig \

bigwigCompare -b1 LNCaP.dht.AR.rep1.-.5end.bigWig \
-b2 LNCaP.dht.AR.rep2.-.5end.bigWig \
--operation add --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed -p 40 \
-bs=1 -o LNCaP.dht.AR.-.5end.bigWig \

bigwigCompare -b1 LNCaP.veh.AR.rep1.+.5end.bigWig \
-b2 LNCaP.veh.AR.rep2.+.5end.bigWig \
--operation add --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed -p 40 \
-bs=1 -o LNCaP.veh.AR.+.5end.bigWig \

bigwigCompare -b1 LNCaP.veh.AR.rep1.-.5end.bigWig \
-b2 LNCaP.veh.AR.rep2.-.5end.bigWig \
--operation add --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed -p 40 \
-bs=1 -o LNCaP.veh.AR.-.5end.bigWig \

bigwigCompare -b1 malignant.1.AR.rep1.+.5end.bigWig \
-b2 malignant.1.AR.rep2.+.5end.bigWig \
--operation add --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed -p 40 \
-bs=1 -o malignant.1.AR.+.5end.bigWig \

bigwigCompare -b1 malignant.1.AR.rep1.-.5end.bigWig \
-b2 malignant.1.AR.rep2.-.5end.bigWig \
--operation add --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed -p 40 \
-bs=1 -o malignant.1.AR.-.5end.bigWig \

bigwigCompare -b1 malignant.2.AR.rep1.+.5end.bigWig \
-b2 malignant.2.AR.rep2.+.5end.bigWig \
--operation add --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed -p 40 \
-bs=1 -o malignant.2.AR.+.5end.bigWig \

bigwigCompare -b1 malignant.2.AR.rep1.-.5end.bigWig \
-b2 malignant.2.AR.rep2.-.5end.bigWig \
--operation add --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed -p 40 \
-bs=1 -o malignant.2.AR.-.5end.bigWig \

bigwigCompare -b1 malignant.3.AR.rep1.+.5end.bigWig \
-b2 malignant.3.AR.rep2.+.5end.bigWig \
--operation add --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed -p 40 \
-bs=1 -o malignant.3.AR.+.5end.bigWig \

bigwigCompare -b1 malignant.3.AR.rep1.-.5end.bigWig \
-b2 malignant.3.AR.rep2.-.5end.bigWig \
--operation add --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed -p 40 \
-bs=1 -o malignant.3.AR.-.5end.bigWig \

bigwigCompare -b1 malignant.4.AR.rep1.+.5end.bigWig \
-b2 malignant.4.AR.rep2.+.5end.bigWig \
--operation add --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed -p 40 \
-bs=1 -o malignant.4.AR.+.5end.bigWig \

bigwigCompare -b1 malignant.4.AR.rep1.-.5end.bigWig \
-b2 malignant.4.AR.rep2.-.5end.bigWig \
--operation add --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed -p 40 \
-bs=1 -o malignant.4.AR.-.5end.bigWig \

bigwigCompare -b1 non-malignant.1.AR.rep1.+.5end.bigWig \
-b2 non-malignant.1.AR.rep2.+.5end.bigWig \
--operation add --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed -p 40 \
-bs=1 -o non-malignant.1.AR.+.5end.bigWig \

bigwigCompare -b1 non-malignant.1.AR.rep1.-.5end.bigWig \
-b2 non-malignant.1.AR.rep2.-.5end.bigWig \
--operation add --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed -p 40 \
-bs=1 -o non-malignant.1.AR.-.5end.bigWig \

bigwigCompare -b1 non-malignant.2.AR.rep1.+.5end.bigWig \
-b2 non-malignant.2.AR.rep2.+.5end.bigWig \
--operation add --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed -p 40 \
-bs=1 -o non-malignant.2.AR.+.5end.bigWig \

bigwigCompare -b1 non-malignant.2.AR.rep1.-.5end.bigWig \
-b2 non-malignant.2.AR.rep2.-.5end.bigWig \
--operation add --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed -p 40 \
-bs=1 -o non-malignant.2.AR.-.5end.bigWig \
