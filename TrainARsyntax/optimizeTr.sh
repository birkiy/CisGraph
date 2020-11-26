


http://mitra.stanford.edu/kundaje/avsec/chipnexus/paper/data/chip-nexus/Oct4/counts.pos.subset.bw
http://mitra.stanford.edu/kundaje/avsec/chipnexus/paper/data/chip-nexus/Oct4/counts.neg.subset.bw

../data/chip-nexus/Oct4/counts.pos.subset.bw


cut -f1,2,3 /groups/lackgrp/ll_members/berkay/TermProjectComp541/ARsyntax/results/peak/idr/22RV1.AR-C19.rep.summits.bed --output-delimiter=":"

cut -f1,2,3 /groups/lackgrp/ll_members/berkay/TermProjectComp541/ARsyntax/results/peak/idr/22RV1.AR-C19.rep.summits.bed --output-delimiter=":" | tr : '\t'  > optimize.summit.bed



awk '{print $1"\t"$2"\t"$3}' /groups/lackgrp/ll_members/berkay/TermProjectComp541/ARsyntax/results/peak/idr/22RV1.AR-C19.rep.summits.bed | head


fasta_file: hg19.fa  # reference genome fasta file
task_specs:

  22RV1.AR-C19:
    tracks:
      - 22RV1.AR-C19.merged.+.5end.bigWig
      - 22RV1.AR-C19.merged.-.5end.bigWig
    peaks: optimize.summit.bed




sed -i '/_/d' ./optimize.summit.bed









##############################3


mkdir ARtrainOptimize
cd ARtrainOptimize

bpnet dataspec-stats dataspec.yaml > dataspecStats.log


#!/bin/bash
#SBATCH --job-name=train
#SBATCH --time=23:00:00
#SBATCH --mem-per-cpu=14G
#SBATCH --cpus-per-task=50
#SBATCH --export=all
#SBATCH -p long

expDir=/groups/lackgrp/ll_members/berkay/TermProjectComp541/TrainARsyntax/ARtrainOptimize


uuid=$(uuidgen)
runId=$uuid

echo "==========================================================================\n\n"
echo $runId
echo "==========================================================================\n\n"
modelDir=$expDir/$runId
contribFile=$modelDir/contrib.deeplift.h5
contribNullFile=$modelDir/contrib.deeplift.null.h5
modiscoDir=$modelDir/modisco


echo "==========================================================================\n\n"
echo "TRAIN STARTS"
echo "==========================================================================\n\n"
bpnet train dataspec2.yaml --vmtouch \
  --override="lambda=8.70" \
  . --run-id=$runId --num-workers 50 --in-memory


##################################################################
uuid=8ad2ab4a-59c4-421f-a6a0-e02d6f890f6d


#!/bin/bash
#SBATCH --job-name=contrib1
#SBATCH --time=23:00:00
#SBATCH --mem-per-cpu=14G
#SBATCH --cpus-per-task=50
#SBATCH --export=all
#SBATCH -p long

expDir=/groups/lackgrp/ll_members/berkay/TermProjectComp541/TrainARsyntax/ARtrainOptimize


uuid=8ad2ab4a-59c4-421f-a6a0-e02d6f890f6d
runId=$uuid

echo "==========================================================================\n\n"
echo $runId
echo "==========================================================================\n\n"
modelDir=$expDir/$runId
contribFile=$modelDir/contrib.deeplift.h5
contribNullFile=$modelDir/contrib.deeplift.null.h5
modiscoDir=$modelDir/modisco


echo "==========================================================================\n\n"
echo "CONTRIBUTION SCORES"
echo "==========================================================================\n\n"
bpnet contrib $modelDir --method=deeplift --contrib-wildcard='*/profile/wn' $contribFile --num-workers 50
echo "==========================================================================\n\n"
echo "NULL CONTRIBUTION SCORES"
echo "==========================================================================\n\n"
bpnet contrib $modelDir --method=deeplift --shuffle-seq --max-regions 5000 --contrib-wildcard='*/profile/wn' $contribNullFile --num-workers 50


###############################################3
uuid=cacb676c-fe56-4736-b38a-00e2e2f8ef93

#!/bin/bash
#SBATCH --job-name=contrib2
#SBATCH --time=23:00:00
#SBATCH --mem-per-cpu=14G
#SBATCH --cpus-per-task=50
#SBATCH --export=all
#SBATCH -p long

expDir=/groups/lackgrp/ll_members/berkay/TermProjectComp541/TrainARsyntax/ARtrainOptimize


uuid=cacb676c-fe56-4736-b38a-00e2e2f8ef93
runId=$uuid

echo "==========================================================================\n\n"
echo $runId
echo "==========================================================================\n\n"
modelDir=$expDir/$runId
contribFile=$modelDir/contrib.deeplift.h5
contribNullFile=$modelDir/contrib.deeplift.null.h5
modiscoDir=$modelDir/modisco


echo "==========================================================================\n\n"
echo "CONTRIBUTION SCORES"
echo "==========================================================================\n\n"
bpnet contrib $modelDir --method=deeplift --contrib-wildcard='*/profile/wn' $contribFile --num-workers 50
echo "==========================================================================\n\n"
echo "NULL CONTRIBUTION SCORES"
echo "==========================================================================\n\n"
bpnet contrib $modelDir --method=deeplift --shuffle-seq --max-regions 5000 --contrib-wildcard='*/profile/wn' $contribNullFile --num-workers 50

###################################################






echo "==========================================================================\n\n"
echo "MODISCO RUN"
echo "==========================================================================\n\n"
echo $task;
bpnet modisco-run $contribFile --null-contrib-file=$contribNullFile \
  --contrib-wildcard=$task/profile/wn \
  --premade=modisco-50k --override='TfModiscoWorkflow.min_metacluster_size=1000' \
  --only-task-regions \
  $modiscoDir/$task \
  --num-workers 35 \
  --overwrite


tasks=(
'22RV1.AR-C19'
'22RV1.AR-V7'
'LN95.AR-C19'
'LN95.AR-V7'
'LNCaP.dht.AR'
'LNCaP.veh.AR'
'malignant.1.AR'
'malignant.2.AR'
'malignant.3.AR'
'malignant.4.AR'
'non-malignant.1.AR'
'non-malignant.2.AR'
)

# Run modisco only for the Nanog task
for task in ${tasks[@]}
do
  echo $task;
  bpnet modisco-run $contribFile --null-contrib-file=$contribNullFile \
    --contrib-wildcard=$task/profile/wn \
    --premade=modisco-50k --override='TfModiscoWorkflow.min_metacluster_size=1000' \
    --only-task-regions \
    $modiscoDir/$task \
    --overwrite
done





#!/bin/bash
#SBATCH --job-name=mod
#SBATCH --time=23:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=35
#SBATCH --export=all
#SBATCH -p long


declare -A dataspecs

dataspecs[22RV1.AR-C19]=33.93
dataspecs[22RV1.AR-V7]=32.90
dataspecs[LN95.AR-C19]=17.50
dataspecs[LN95.AR-V7]=25.17
dataspecs[LNCaP.dht.AR]=29.46
dataspecs[LNCaP.veh.AR]=24.68
dataspecs[malignant.1.AR]=13.64
dataspecs[malignant.2.AR]=12.38
dataspecs[malignant.3.AR]=16.48
dataspecs[malignant.4.AR]=12.69
dataspecs[non-malignant.1.AR]=14.66
dataspecs[non-malignant.2.AR]=12.65


expDir=/groups/lackgrp/ll_members/berkay/TermProjectComp541/TrainARsyntax/ARtrainOptimize/dataspecs

for task in ${!dataspecs[@]}
do
  modelDir=$expDir/$task
  contribFile=$modelDir/contrib.deeplift.h5
  contribNullFile=$modelDir/contrib.deeplift.null.h5
  modiscoDir=$modelDir/modisco
  echo "==========================================================================\n\n"
  echo "TRAIN STARTS for $task with lambda ${dataspecs[$task]}"
  echo "==========================================================================\n\n"
  bpnet train $task".yaml" --vmtouch \
    --override="lambda=${dataspecs[$task]}" \
    . --run-id=$task --num-workers 35 --in-memory


  echo "==========================================================================\n\n"
  echo "CONTRIBUTION SCORES"
  echo "==========================================================================\n\n"
  bpnet contrib $modelDir --method=deeplift --contrib-wildcard='*/profile/wn' $contribFile --num-workers 35
  echo "==========================================================================\n\n"
  echo "NULL CONTRIBUTION SCORES"
  echo "==========================================================================\n\n"
  bpnet contrib $modelDir --method=deeplift --shuffle-seq --max-regions 5000 --contrib-wildcard='*/profile/wn' $contribNullFile --num-workers 35


  echo "==========================================================================\n\n"
  echo "MODISCO RUN"
  echo "==========================================================================\n\n"
  echo $task;
  bpnet modisco-run $contribFile --null-contrib-file=$contribNullFile \
    --contrib-wildcard=$task/profile/wn \
    --premade=modisco-50k --override='TfModiscoWorkflow.min_metacluster_size=1000' \
    --only-task-regions \
    $modiscoDir/$task \
    --num-workers 35 \
    --overwrite
done


bpnet modisco-run contrib.deeplift.h5 --null-contrib-file=contrib.deeplift.null.h5 \
  --contrib-wildcard=22RV1.AR-C19/profile/wn \
  --premade=modisco-50k --override='TfModiscoWorkflow.min_metacluster_size=1000' \
  --only-task-regions \
  modisco/22RV1.AR-C19 \
  --num-workers 35 \
  --overwrite





import seaborn as sns
from bpnet.cli.contrib import ContribFile
from bpnet.plot.tracks import plot_tracks, to_neg
import seaborn as sns
import matplotlib.pyplot as plt

cf = ContribFile("8eb61138-ad96-4c71-bb7f-b641baa6b93d/deeplift")

profiles = cf.get_profiles()
contrib_scores = cf.get_contrib()

examples = list({v.max(axis=-2).mean(axis=-1).argmax() for k,v in profiles.items()})
examples

from bpnet.dataspecs import DataSpec, TaskSpec
from pathlib import Path


model_dir = "/groups/lackgrp/ll_members/berkay/TermProjectComp541/TrainARsyntax/ARtrainOptimize/8eb61138-ad96-4c71-bb7f-b641baa6b93d"
model_dir = Path(model_dir)
ds = DataSpec.load(model_dir / 'dataspec.yml')
tasks = list(ds.task_specs)



xrange = slice(50, 150)
for idx in examples:
  plot_tracks({**{'profile/' + k: to_neg(v[idx,xrange]) for k,v in profiles.items()},
               **{'contrib/' + k:v[idx,xrange] for k,v in contrib_scores.items()}},
             title=idx,
             rotate_y=0,
             fig_width=10,
             fig_height_per_track=1);
  sns.despine(top=True, right=True, bottom=True)


fig.savefig("contrib.pdf")



bpnet export-bw . --regions=intervals.bed --scale-contribution bigwigs/

bpnet modisco-run contrib.scores.h5 --premade=modisco-50k --override='TfModiscoWorkflow.max_seqlets_per_metacluster=20000' modisco/

bpnet cwm-scan modisco/ --contrib-file=contrib.scores.h5 modisco/motif-instances.tsv.gz
