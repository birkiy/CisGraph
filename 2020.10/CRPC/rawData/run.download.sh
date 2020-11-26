#!/bin/bash
#SBATCH --job-name=crpc.download
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=30
#SBATCH --export=all
#SBATCH -p long

echo "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48308"

parallel-fastq-dump -t 30 --gzip --split-3 -s SRR922203
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR922204
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR922206
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR922207
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR922209
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR922210
parallel-fastq-dump -t 30 --gzip --split-3 -s SRR922212
