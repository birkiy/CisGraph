

colorPalette = {"cre": "#F5B278",
                "tad": "#40A59A",
                "com": "#EA846A",
                "chr": "#395560"}


server = True

if server:
    home = "/kuacc/users/ualtintas20"
    envR = "/kuacc/apps/R/3.6.1/bin/Rscript"
else:
    home = "/home/birkiy"
    envR = f"{home}/anaconda3/envs/CisGraph/bin/Rscript"

projectRoot = f"{home}/github/CisGraph/Vers1.0"
dataRoot = f"{home}/github/Data/CisGraph/Vers1.0"
figureRoot = f"{home}/github/Figures/CisGraph/Vers1.0"



import csv
import pickle
import subprocess
import os

import numpy as np
import networkx as nx
import pandas as pd
from scipy import sparse as sps

from matplotlib import pyplot, patches
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import gridspec, colors
# from pygraphviz import *
from matplotlib.colors import LinearSegmentedColormap
from statannot import add_stat_annotation

import matplotlib
import itertools
import scipy
import scipy.stats
import math
import random
import heapq
from collections import deque
from statsmodels.sandbox.stats.multicomp import multipletests
from statannot import add_stat_annotation

from _collections import deque
import collections
import bisect as bi
import time
