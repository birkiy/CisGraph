
import bpnet
from bpnet.cli.contrib import ContribFile
from bpnet.plot.tracks import plot_tracks, to_neg

import uuid
import datetime
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import subprocess


import os


from pathlib import Path
exp_dir = Path("/groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/train")



os.system(f"head {exp_dir}/dataspec.yml")

os.system(f"cat {exp_dir}/config.gin")




model_dir = exp_dir / 'output'
contrib_file = model_dir/'contrib.deeplift.h5'
contrib_null_file = model_dir/'contrib.deeplift.null.h5'
modisco_dir = model_dir/'modisco'


run_id = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S") + "_" + str(uuid.uuid4())
run_id = run_id + ".exampleRun"

os.system(f"cd {exp_dir} && bpnet train dataspec.yml --premade=bpnet9 --config=config.gin . --override='train.epochs=10' --run-id '{run_id}' --in-memory")


os.system(f"rm {exp_dir}/output && ln -srf {exp_dir}/{run_id}  {exp_dir}/output")


HTML(filename=f"{exp_dir}/{run_id}/evaluate.html")

# subprocess.run(f"cat {exp_dir}/dataspec.yml")
#
# tasks = ["AR.1", "FOXA1.1"]
#
# for task in tasks:
#     print(task)
#
#     # shell("zcat bpnet/examples/data/chip-nexus/{task}/idr-optimal-set.summit.subset.bed.gz | wc -l")7




base_url = 'http://mitra.stanford.edu/kundaje/avsec/chipnexus/paper/data'

dataspecs = ['chip-nexus/dataspec.yml',
             'chip-seq/dataspec.yml']

paths = ["counts.neg.subset.bw", "counts.pos.subset.bw",  "idr-optimal-set.summit.subset.bed.gz"]
tfs = ["Nanog", "Oct4", "Sox2"]

fRun = [f"wget {base_url}/{tf}/{path} -O data/{tf}/{path}" for tf in tfs for path in paths]

paths = ["counts.neg.subset.bw", "counts.pos.subset.bw"]
fRun += [f"wget {base_url}/patchcap/{path} -O data/patchcap/{path}" for path in paths]


for r in fRun:
    os.system(r)




from pathlib import Path
exp_dir = Path('bpnet/examples/chip-nexus/')




os.system(f"cat {exp_dir}/dataspec.yml")

os.system(f"zcat bpnet/examples/data/chip-nexus/*/idr-optimal-set.summit.subset.bed.gz | cut -f 1 | sort -u")

tasks = ['Oct4', 'Sox2', 'Nanog']


for task in tasks:
    print(task)
    os.system(f"zcat bpnet/examples/data/chip-nexus/{task}/idr-optimal-set.summit.subset.bed.gz | wc -l")


os.system(f"cat {exp_dir}/config.gin")

model_dir = exp_dir / 'output'
contrib_file = model_dir/'contrib.deeplift.h5'
contrib_null_file = model_dir/'contrib.deeplift.null.h5'
modisco_dir = model_dir/'modisco'


# setup a new run_id (could be done automatically, but then the output directory would change)
run_id = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S") + "_" + str(uuid.uuid4())

# Train for at most 10 epochs
os.system(f"cd {exp_dir} && bpnet train dataspec.yml --premade=bpnet9 --config=config.gin . --override='train.epochs=10' --run-id '{run_id}' --in-memory")

# softlink the new output directory
os.system(f"rm {exp_dir}/output && ln -srf {exp_dir}/{run_id}  {exp_dir}/output")

model_dir = f"{exp_dir}/{run_id}"


run_id = "2020-10-03_07-55-17_1e61b98d-8bb0-4220-b501-3def6877fd00"

# contribution scores
os.system(f"bpnet contrib {model_dir} --method=deeplift --memfrac-gpu=1 --contrib-wildcard='*/profile/wn' {contrib_file}")




import seaborn as sns
from bpnet.cli.contrib import ContribFile
from bpnet.plot.tracks import plot_tracks, to_neg
import seaborn as sns
import matplotlib.pyplot as plt

cf = ContribFile(contrib_file)

profiles = cf.get_profiles()
contrib_scores = cf.get_contrib()

examples = list({v.max(axis=-2).mean(axis=-1).argmax() for k,v in profiles.items()})
examples

tasks = ['Oct4', 'Sox2', 'Nanog']



fig = plt.figure(figsize=[8,8])

xrange = slice(50, 150)
for idx in examples:
    plot_tracks({**{'profile/' + k: to_neg(v[idx,xrange]) for k,v in profiles.items()},
                **{'contrib/' + k:v[idx,xrange] for k,v in contrib_scores.items()}},
                title=idx,
                rotate_y=0,
                fig_width=10,
                fig_height_per_track=1);
    sns.despine(top=True, right=True, bottom=True)


fig.savefig(f"/groups/lackgrp/ll_members/berkay/SYNTAX/exampleRun/results/train/bpnet/examples/chip-nexus/output/contributions.pdf")
