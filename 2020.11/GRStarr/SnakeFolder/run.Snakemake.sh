#!/bin/bash


snakemake \
  --snakefile code/Snakefile \
  --configfile code/config.yaml \
  --unlock -j1


# snakemake version: 5.26.1
# conda activate StarrAnalysis

snakemake \
  --snakefile code/Snakefile \
  --configfile code/config.yaml \
  -j100 \
  --cluster-config code/cluster.yaml  \
  --rerun-incomplete \
  --use-conda \
  --cluster "sbatch \
      -A {cluster.partition} \
      -c {cluster.c} \
      -t {cluster.time} \
      --mem {cluster.mem} \
      -o logsSlurm/{rule}_{wildcards} \
      -e logsSlurm/{rule}_{wildcards} "



snakemake \
  --snakefile code/Snakefile \
  --configfile code/config.yaml \
  -n --debug-dag -r


snakemake \
  --snakefile code/Snakefile \
  --configfile code/config.yaml \
  -n --dag | dot -Tpdf > StarrSamples.pdf


snakemake \
  --snakefile Snakefile \
  --configfile config.yaml \
  -n --forceall --rulegraph | dot -Tpdf > ARsyntaxRule.pdf
