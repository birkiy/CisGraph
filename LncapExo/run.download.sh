#!/bin/bash
#SBATCH --job-name=star.pbs
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=10
#SBATCH --export=all
#SBATCH -p long



parallel-fastq-dump -t 10 --gzip -O SRX5658368	-s SRR8871335
parallel-fastq-dump -t 10 --gzip -O SRX5658367	-s SRR8871336
parallel-fastq-dump -t 10 --gzip -O SRX5658366	-s SRR8871337
parallel-fastq-dump -t 10 --gzip -O SRX5658365	-s SRR8871338
parallel-fastq-dump -t 10 --gzip -O SRX5658364	-s SRR8871339
parallel-fastq-dump -t 10 --gzip -O SRX5658363	-s SRR8871340
parallel-fastq-dump -t 10 --gzip -O SRX5658362	-s SRR8871341
parallel-fastq-dump -t 10 --gzip -O SRX5658361	-s SRR8871342
parallel-fastq-dump -t 10 --gzip -O SRX5658360	-s SRR8871343
parallel-fastq-dump -t 10 --gzip -O SRX5658359	-s SRR8871344
parallel-fastq-dump -t 10 --gzip -O SRX5658358	-s SRR8871345
parallel-fastq-dump -t 10 --gzip -O SRX5658357	-s SRR8871346
parallel-fastq-dump -t 10 --gzip -O SRX5658356	-s SRR8871347
parallel-fastq-dump -t 10 --gzip -O SRX5658355	-s SRR8871348
parallel-fastq-dump -t 10 --gzip -O SRX5658354	-s SRR8871349
parallel-fastq-dump -t 10 --gzip -O SRX5658353	-s SRR8871350
parallel-fastq-dump -t 10 --gzip -O SRX5658352	-s SRR8871351
parallel-fastq-dump -t 10 --gzip -O SRX5658351	-s SRR8871352
parallel-fastq-dump -t 10 --gzip -O SRX5658350	-s SRR8871353
parallel-fastq-dump -t 10 --gzip -O SRX5658349	-s SRR8871354
parallel-fastq-dump -t 10 --gzip -O SRX5658348	-s SRR8871355
parallel-fastq-dump -t 10 --gzip -O SRX5658347	-s SRR8871356
parallel-fastq-dump -t 10 --gzip -O SRX5658346	-s SRR8871357
parallel-fastq-dump -t 10 --gzip -O SRX5658345	-s SRR8871358
parallel-fastq-dump -t 10 --gzip -O SRX5658344	-s SRR8871359
parallel-fastq-dump -t 10 --gzip -O SRX5658343	-s SRR8871360
parallel-fastq-dump -t 10 --gzip -O SRX5658342	-s SRR8871361
parallel-fastq-dump -t 10 --gzip -O SRX5658341	-s SRR8871362
parallel-fastq-dump -t 10 --gzip -O SRX5658340	-s SRR8871363
parallel-fastq-dump -t 10 --gzip -O SRX5658339	-s SRR8871364
parallel-fastq-dump -t 10 --gzip -O SRX5658338	-s SRR8871365
parallel-fastq-dump -t 10 --gzip -O SRX5658337	-s SRR8871366
parallel-fastq-dump -t 10 --gzip -O SRX5658336	-s SRR8871367
parallel-fastq-dump -t 10 --gzip -O SRX5658335	-s SRR8871368
parallel-fastq-dump -t 10 --gzip -O SRX5658334	-s SRR8871369
parallel-fastq-dump -t 10 --gzip -O SRX5658333	-s SRR8871370





# parallel-fastq-dump -t 10 --gzip -O SRX882919	-s SRR1810086
# parallel-fastq-dump -t 10 --gzip -O SRX882920	-s SRR1810087
# parallel-fastq-dump -t 10 --gzip -O SRX882921	-s SRR1810088
# parallel-fastq-dump -t 10 --gzip -O SRX882922	-s SRR1810089
# parallel-fastq-dump -t 10 --gzip -O SRX882923	-s SRR1810090
# parallel-fastq-dump -t 10 --gzip -O SRX882924	-s SRR1810091
# parallel-fastq-dump -t 10 --gzip -O SRX882925	-s SRR1810092
# parallel-fastq-dump -t 10 --gzip -O SRX882926	-s SRR1810093
# parallel-fastq-dump -t 10 --gzip -O SRX882927	-s SRR1810094
# parallel-fastq-dump -t 10 --gzip -O SRX882928	-s SRR1810095
# parallel-fastq-dump -t 10 --gzip -O SRX882929	-s SRR1810096
# parallel-fastq-dump -t 10 --gzip -O SRX882930	-s SRR1810097
# parallel-fastq-dump -t 10 --gzip -O SRX882931	-s SRR1810098
# parallel-fastq-dump -t 10 --gzip -O SRX882932	-s SRR1810099
# parallel-fastq-dump -t 10 --gzip -O SRX882933	-s SRR1810100
# parallel-fastq-dump -t 10 --gzip -O SRX882934	-s SRR1810101
# parallel-fastq-dump -t 10 --gzip -O SRX882935	-s SRR1810102
# parallel-fastq-dump -t 10 --gzip -O SRX882936	-s SRR1810103
# parallel-fastq-dump -t 10 --gzip -O SRX882937	-s SRR1810104
# parallel-fastq-dump -t 10 --gzip -O SRX882938	-s SRR1810105
# parallel-fastq-dump -t 10 --gzip -O SRX882939	-s SRR1810106
# parallel-fastq-dump -t 10 --gzip -O SRX882940	-s SRR1810107
# parallel-fastq-dump -t 10 --gzip -O SRX882941	-s SRR1810108




parallel-fastq-dump -t 20 --gzip -O SRX882911 -s SRR1810079
parallel-fastq-dump -t 10 --gzip -O SRX882912 -s SRR1810080
parallel-fastq-dump -t 10 --gzip -O SRX882913 -s SRR1810081
