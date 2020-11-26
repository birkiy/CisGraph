#!/bin/bash

sacct --format=User,Account,JobID,Jobname,partition,state,time,start,end,elapsed,MaxRss,MaxVMSize,nnodes,ncpus,nodelis

snakemake \
  --snakefile code/Snakefile \
  --configfile code/config.yaml \
  -j100 \
  --cluster-config code/cluster.yaml \
  --rerun-incomplete \
  --use-conda \
  --cluster "sbatch \
    -A {cluster.account} \
    -p {cluster.partition} \
    -t {cluster.time} \
    --mem {cluster.mem} \
    -c {cluster.c} \
    -o logsSlurm/{rule}_{wildcards} \
    -e logsSlurm/{rule}_{wildcards} \
    --exclude cn13"


sbatch \
  -p normal,big-mem,long \
  -t {resources.time_min} \
  --mem={resources.mem_mb} \
  -c {resources.cpus} \
  -o logs_slurm/{rule}_{wildcards} \
  -e logs_slurm/{rule}_{wildcards}



snakemake \
  --snakefile code/Snakefile \
  --configfile code/config.yaml \
  --profile slurm





{
  "__default__" :
  {
      "account" : "ualtintas",
      "time" : "23:59:00",
      "partition" :"normal,big-mem,long"
  },
  "jobs" : "100",
  "default-resources" : { "cpus": "1", "mem" : "1G"}
}



account: ualtintas
time: 23:59:00
partition:
  - normal
  - big-mem
  - long
jobs: 100
cluster: sbatch \
  -A ualtintas \
  -p normal, big-mem, long \
  -t 23:59:00 \
  --mem={resources.mem} \
  -c {resources.cpus} \
  -o logs_slurm/{rule}_{wildcards} \
  -e logs_slurm/{rule}_{wildcards}
use-conda: true
default-resources:
  - cpus=1
  - mem=1000
resources:
  - cpusHigh=16
  - memHigh=320000
  - cpusHigh=32
  - memHigh=640000



snakemake \
  --snakefile code/Snakefile \
  --configfile code/config.yaml \
  -n --dag | dot -Tpdf > Samples.pdf
