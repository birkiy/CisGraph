#!/bin/bash


snakemake \
  --snakefile code/Snakefile \
  --configfile code/config.yaml \
  --unlock -j1


# snakemake version: 5.26.1

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
      --output 'slurms/%j.out' \
      --error 'slurms/%j.err'"


snakemake \
  --snakefile code/Snakefile \
  --configfile code/config.yaml \
  -n --lt


snakemake \
  --snakefile code/Snakefile \
  --configfile code/config.yaml \
  -n --dag | dot -Tpdf > StarrSamples.pdf


snakemake \
  --snakefile Snakefile \
  --configfile config.yaml \
  -n --forceall --rulegraph | dot -Tpdf > ARsyntaxRule.pdf













 bamCompare -b1 results/mapping/processed/GR.0h.dex.merged.final.bam -b2 results/mapping/processed/input.GSE114063.merged.final.bam -o results/bigwig/GR.0h.dex.merged.RPKM.SES.bigWig                 --blackListFileName /home/ualtintas/genomeAnnotations/ENCFF001TDO.bed                 --extendReads 150                 --centerReads                 --scaleFactorsMethod SES                 --operation ratio                 -p 16                --scaleFactors 1.2:0.5                 --normalizeUsing {wildcards.type}
