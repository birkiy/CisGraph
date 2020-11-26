#!/bin/bash

snakemake \
  --snakefile Snakefile \
  --configfile config.yaml \
  -j100 \
  --cluster-config cluster.yaml  \
  --rerun-incomplete \
  --use-conda \
  --cluster "sbatch \
      -A {cluster.partition} \
      -c {cluster.c} \
      -t {cluster.time} \
      --mem {cluster.mem} \
      --output 'slurms/%j.out' \
      --error 'slurms/%j.err'"



snakemake \
  --snakefile Snakefile \
  --configfile config.yaml \
  -n --dag | dot -Tpdf > CRPCSamples.pdf


snakemake \
  --snakefile Snakefile \
  --configfile config.yaml \
  --unlock

snakemake \
  --snakefile Snakefile \
  --configfile config.yaml \
  -n --forceall --rulegraph | dot -Tpdf > CRPCRule.pdf
