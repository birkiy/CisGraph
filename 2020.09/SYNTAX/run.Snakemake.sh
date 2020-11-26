#!/bin/bash
#SBATCH --job-name=Chip.download
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=64
#SBATCH --export=all
#SBATCH -p long

snakemake --cores 64
