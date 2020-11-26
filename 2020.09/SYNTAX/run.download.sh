#!/bin/bash
#SBATCH --job-name=Chip.download
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=30
#SBATCH --export=all
#SBATCH -p long

echo "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94682"

parallel-fastq-dump -t 30 --gzip --split-3 -s SRR5238072
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR5238073
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR5238074
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR5238075
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR5238076
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR5238077
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR5238079
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR5238080
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR5238081
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR5238084
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR5238086
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR5238089
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR5238090
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR5238093
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR5238095
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR5238096
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR5238097
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR5238098
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR5238101
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR5238103
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR5238104
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR5238105
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR5238083
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR5238078
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR5238088
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR5238091
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR5238100
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR5238082
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR5238085
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR5238087
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR5238092
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR5238099
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR5238102
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR5238094



mv SRR5238072.fastq.gz rawData/AR.veh.rep1.fastq.gz
mv SRR5238073.fastq.gz rawData/AR.r1881.rep1.fastq.gz
mv SRR5238074.fastq.gz rawData/AR.veh.rep2.fastq.gz
mv SRR5238075.fastq.gz rawData/AR.r1881.rep2.fastq.gz
mv SRR5238090.fastq.gz rawData/FOXA1.veh.rep1.fastq.gz
mv SRR5238093.fastq.gz rawData/FOXA1.r1881.rep1.fastq.gz
mv SRR5238095.fastq.gz rawData/FOXA1.veh.rep2.fastq.gz
mv SRR5238096.fastq.gz rawData/FOXA1.r1881.rep2.fastq.gz



mv SRR5238076.fastq.gz rawData/ARID.veh.rep1.fastq.gz
mv SRR5238077.fastq.gz rawData/ARID.r1881.rep1.fastq.gz
mv SRR5238079.fastq.gz rawData/ARID.veh.rep2.fastq.gz
mv SRR5238080.fastq.gz rawData/ARID.r1881.rep2.fastq.gz
mv SRR5238081.fastq.gz rawData/BRG1.veh.rep1.fastq.gz
mv SRR5238084.fastq.gz rawData/BRG1.r1881.rep1.fastq.gz
mv SRR5238086.fastq.gz rawData/BRG1.veh.rep2.fastq.gz
mv SRR5238089.fastq.gz rawData/BRG1.r1881.rep2.fastq.gz
mv SRR5238097.fastq.gz rawData/HOB13.veh.rep1.fastq.gz
mv SRR5238098.fastq.gz rawData/HOB13.r1881.rep1.fastq.gz
mv SRR5238101.fastq.gz rawData/HOB13.veh.rep2.fastq.gz
mv SRR5238103.fastq.gz rawData/HOB13.r1881.rep2.fastq.gz
mv SRR5238104.fastq.gz rawData/HOB13.veh.rep3.fastq.gz
mv SRR5238105.fastq.gz rawData/HOB13.r1881.rep3.fastq.gz
mv SRR5238083.fastq.gz rawData/TLE3.veh.rep1.fastq.gz
mv SRR5238078.fastq.gz rawData/TLE3.r1881.rep1.fastq.gz
mv SRR5238088.fastq.gz rawData/TLE3.veh.rep2.fastq.gz
mv SRR5238091.fastq.gz rawData/TLE3.r1881.rep2.fastq.gz
mv SRR5238100.fastq.gz rawData/TRIM28.veh.rep1.fastq.gz
mv SRR5238082.fastq.gz rawData/TRIM28.r1881.rep1.fastq.gz
mv SRR5238085.fastq.gz rawData/TRIM28.veh.rep2.fastq.gz
mv SRR5238087.fastq.gz rawData/TRIM28.r1881.rep2.fastq.gz
mv SRR5238092.fastq.gz rawData/WDHD1.veh.rep1.fastq.gz
mv SRR5238099.fastq.gz rawData/WDHD1.r1881.rep1.fastq.gz
mv SRR5238102.fastq.gz rawData/WDHD1.veh.rep2.fastq.gz
mv SRR5238094.fastq.gz rawData/WDHD1.r1881.rep2.fastq.gz



































#!/bin/bash
#SBATCH --job-name=Chip.download
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=30
#SBATCH --export=all
#SBATCH -p long


echo "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE83860"

parallel-fastq-dump -t 30 --gzip --split-3 -s SRR3728821
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR3728822
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR3728823
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR3728824
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR3728825
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR3728826
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR3728828
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR3728827
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR3728837
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR3728838
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR3728839
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR3728840
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR3728841
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR3728842
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR3728843
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR3728844
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR3728865
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR3728866
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR3728867
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR3728868


cat SRR3728821.fastq.gz SRR3728822.fastq.gz > rawData/AR.dht.rep1.fastq.gz
cat SRR3728823.fastq.gz SRR3728824.fastq.gz > rawData/AR.dht.rep2.fastq.gz
cat SRR3728825.fastq.gz SRR3728826.fastq.gz > rawData/AR.dmso.rep1.fastq.gz
cat SRR3728828.fastq.gz SRR3728827.fastq.gz > rawData/AR.dmso.rep2.fastq.gz
cat SRR3728837.fastq.gz SRR3728838.fastq.gz > rawData/FOXA1.dht.rep1.fastq.gz
cat SRR3728839.fastq.gz SRR3728840.fastq.gz > rawData/FOXA1.dht.rep2.fastq.gz
cat SRR3728841.fastq.gz SRR3728842.fastq.gz > rawData/FOXA1.dmso.rep1.fastq.gz
cat SRR3728843.fastq.gz SRR3728844.fastq.gz > rawData/FOXA1.dmso.rep2.fastq.gz
cat SRR3728865.fastq.gz SRR3728866.fastq.gz > rawData/input.rep1.fastq.gz
cat SRR3728867.fastq.gz SRR3728868.fastq.gz > rawData/input.rep2.fastq.gz


mv rawData/input.rep1.fastq.gz rawData/input.GSE83860.rep1.fastq.gz 
mv rawData/input.rep2.fastq.gz rawData/input.GSE83860.rep2.fastq.gz

#!/bin/bash
#SBATCH --job-name=Chip.download
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=30
#SBATCH --export=all
#SBATCH -p long

echo "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117304"

parallel-fastq-dump -t 30 --gzip --split-3 -s SRR7536832
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR7536834
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR8816206
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR8816207


mv SRR7536832.fastq.gz rawData/input.GSE117304.dht.fastq.gz
mv SRR7536834.fastq.gz rawData/input.GSE117304.etoh.fastq.gz
mv SRR8816206.fastq.gz rawData/HOXB13.dht.rep1.fastq.gz
mv SRR8816207.fastq.gz rawData/HOXB13.etoh.rep1.fastq.gz
