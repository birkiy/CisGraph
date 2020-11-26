import bpnet
from bpnet.cli.contrib import ContribFile
from bpnet.plot.tracks import plot_tracks, to_neg

import uuid
import datetime
import seaborn as sns
import matplotlib.pyplot as plt
from IPython.display import HTML
import pandas as pd
import numpy as np


from pathlib import Path
import os


exDir = Path('/groups/lackgrp/ll_members/berkay/ARsyntax/results/train')


exDir = Path('/groups/lackgrp/ll_members/berkay/ARsyntax/optimize/train')

modelDir = exDir / 'output'
contribFile = modelDir/'contrib.deeplift.h5'
contribNullFile = modelDir/'contrib.deeplift.null.h5'
modiscoDir = modelDir/'modisco'

runId = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S") + "_" + str(uuid.uuid4())
os.system("cat config.gin")


cd {exDir} &&
os.system(f"bpnet train dataspec.yaml --premade=bpnet9 --config='config.gin' . --override='train.epochs=10' --run-id '{runId}' --in-memory")

os.system(f"bpnet train dataspec.yaml --premade=bpnet9 --config='config.gin' . --run-id '{runId}' --in-memory")

bpnet train dataspec.yml --config=config.gin . --override='train.epochs=10' --run-id "firstRun"  --in-memory

bpnet train dataspec.yml .

os.system(f"rm {exDir}/output && ln -srf {exDir}/{runId}  {exDir}/output")









bpnet train dataspec.yaml --premade=bpnet9 --config='config3L50E200S.gin' . --run-id 'con3L50E200S' --in-memory"
