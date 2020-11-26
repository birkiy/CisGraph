#!/bin/bash
#SBATCH --job-name=star.download
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=30
#SBATCH --export=all
#SBATCH -p long

# parallel-fastq-dump -t 30 --gzip -s SRR7122167 --split-3
# parallel-fastq-dump -t 30 --gzip -s SRR7122172 --split-3
# parallel-fastq-dump -t 30 --gzip -s SRR7122177 --split-3
# parallel-fastq-dump -t 30 --gzip -s SRR7122182 --split-3
# parallel-fastq-dump -t 30 --gzip -s SRR7122187 --split-3
# parallel-fastq-dump -t 30 --gzip -s SRR7122192 --split-3
# parallel-fastq-dump -t 30 --gzip -s SRR7122193 --split-3
# parallel-fastq-dump -t 30 --gzip -s SRR7122194 --split-3
# parallel-fastq-dump -t 30 --gzip -s SRR7122195 --split-3
# parallel-fastq-dump -t 30 --gzip -s SRR7122197 --split-3
# parallel-fastq-dump -t 30 --gzip -s SRR7122198 --split-3
# parallel-fastq-dump -t 30 --gzip -s SRR7122196 --split-3
# parallel-fastq-dump -t 30 --gzip -s SRR7122199 --split-3
# parallel-fastq-dump -t 30 --gzip -s SRR7122200 --split-3
# parallel-fastq-dump -t 30 --gzip -s SRR7122201 --split-3
# parallel-fastq-dump -t 30 --gzip -s SRR7122202 --split-3
# parallel-fastq-dump -t 30 --gzip -s SRR7122203 --split-3



# parallel-fastq-dump -t 30 --gzip --split-3 -s SRR7122168
# parallel-fastq-dump -t 30 --gzip --split-3 -s SRR7122169
# parallel-fastq-dump -t 30 --gzip --split-3 -s SRR7122170
# parallel-fastq-dump -t 30 --gzip --split-3 -s SRR7122171
# parallel-fastq-dump -t 30 --gzip --split-3 -s SRR7122173
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR7122174
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR7122175
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR7122176
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR7122178
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR7122179
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR7122180
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR7122181
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR7122183
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR7122184
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR7122185
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR7122186
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR7122188
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR7122189
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR7122190
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR7122191





# mv SRR7122167_1.fastq.gz	GR.0h.dex.rep1.R1.fastq.gz
# mv SRR7122167_2.fastq.gz	GR.0h.dex.rep1.R2.fastq.gz
# mv SRR7122168_1.fastq.gz	GR.0h.dex.rep2.R1.fastq.gz
# mv SRR7122168_2.fastq.gz	GR.0h.dex.rep2.R2.fastq.gz
# mv SRR7122169_1.fastq.gz	GR.0h.dex.rep3.R1.fastq.gz
# mv SRR7122169_2.fastq.gz	GR.0h.dex.rep3.R2.fastq.gz
# mv SRR7122170_1.fastq.gz	GR.0h.dex.rep4.R1.fastq.gz
# mv SRR7122170_2.fastq.gz	GR.0h.dex.rep4.R2.fastq.gz
# mv SRR7122171_1.fastq.gz	GR.0h.dex.rep5.R1.fastq.gz
# mv SRR7122171_2.fastq.gz	GR.0h.dex.rep5.R2.fastq.gz
# mv SRR7122172_1.fastq.gz	GR.1h.dex.rep1.R1.fastq.gz
# mv SRR7122172_2.fastq.gz	GR.1h.dex.rep1.R2.fastq.gz
# mv SRR7122173_1.fastq.gz	GR.1h.dex.rep2.R1.fastq.gz
# mv SRR7122173_2.fastq.gz	GR.1h.dex.rep2.R2.fastq.gz
# mv SRR7122174_1.fastq.gz	GR.1h.dex.rep3.R1.fastq.gz
# mv SRR7122174_2.fastq.gz	GR.1h.dex.rep3.R2.fastq.gz
# mv SRR7122175_1.fastq.gz	GR.1h.dex.rep4.R1.fastq.gz
# mv SRR7122175_2.fastq.gz	GR.1h.dex.rep4.R2.fastq.gz
# mv SRR7122177_1.fastq.gz	GR.4h.dex.rep1.R1.fastq.gz
# mv SRR7122177_2.fastq.gz	GR.4h.dex.rep1.R2.fastq.gz
# mv SRR7122182_1.fastq.gz	GR.8h.dex.rep1.R1.fastq.gz
# mv SRR7122182_2.fastq.gz	GR.8h.dex.rep1.R2.fastq.gz
# mv SRR7122187_1.fastq.gz	GR.12h.dex.rep1.R1.fastq.gz
# mv SRR7122187_2.fastq.gz	GR.12h.dex.rep1.R2.fastq.gz
# mv SRR7122176_1.fastq.gz GR.1h.dex.rep5.R1.fastq.gz
# mv SRR7122176_2.fastq.gz GR.1h.dex.rep5.R2.fastq.gz



# mv SRR7122178_1.fastq.gz	GR.4h.dex.rep2.R1.fastq.gz
# mv SRR7122178_2.fastq.gz	GR.4h.dex.rep2.R2.fastq.gz
# mv SRR7122179_1.fastq.gz	GR.4h.dex.rep3.R1.fastq.gz
# mv SRR7122179_2.fastq.gz	GR.4h.dex.rep3.R2.fastq.gz
# mv SRR7122180_1.fastq.gz	GR.4h.dex.rep4.R1.fastq.gz
# mv SRR7122180_2.fastq.gz	GR.4h.dex.rep4.R2.fastq.gz
# mv SRR7122181_1.fastq.gz	GR.4h.dex.rep5.R1.fastq.gz
# mv SRR7122181_2.fastq.gz	GR.4h.dex.rep5.R2.fastq.gz
# mv SRR7122183_1.fastq.gz	GR.8h.dex.rep2.R1.fastq.gz
# mv SRR7122183_2.fastq.gz	GR.8h.dex.rep2.R2.fastq.gz
# mv SRR7122184_1.fastq.gz	GR.8h.dex.rep3.R1.fastq.gz
# mv SRR7122184_2.fastq.gz	GR.8h.dex.rep3.R2.fastq.gz
# mv SRR7122185_1.fastq.gz	GR.8h.dex.rep4.R1.fastq.gz
# mv SRR7122185_2.fastq.gz	GR.8h.dex.rep4.R2.fastq.gz
# mv SRR7122186_1.fastq.gz	GR.8h.dex.rep5.R1.fastq.gz
# mv SRR7122186_2.fastq.gz	GR.8h.dex.rep5.R2.fastq.gz
# mv SRR7122188_1.fastq.gz	GR.12h.dex.rep2.R1.fastq.gz
# mv SRR7122188_2.fastq.gz	GR.12h.dex.rep2.R2.fastq.gz
# mv SRR7122189_1.fastq.gz	GR.12h.dex.rep3.R1.fastq.gz
# mv SRR7122189_2.fastq.gz	GR.12h.dex.rep3.R2.fastq.gz
# mv SRR7122190_1.fastq.gz	GR.12h.dex.rep4.R1.fastq.gz
# mv SRR7122190_2.fastq.gz	GR.12h.dex.rep4.R2.fastq.gz
# mv SRR7122191_1.fastq.gz	GR.12h.dex.rep5.R1.fastq.gz
# mv SRR7122191_2.fastq.gz	GR.12h.dex.rep5.R2.fastq.gz
