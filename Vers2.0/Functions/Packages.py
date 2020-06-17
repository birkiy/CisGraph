

colorPalette = {"cre": "#F5B278",
                "tad": "#40A59A",
                "com": "#EA846A",
                "chr": "#395560"}


server = True

if server:
    home = "/home/ualtintas/github/CisGraph/Vers2.0"
    envR = "/home/ualtintas/anaconda3/envs/CisGraph/bin/Rscript"
else:
    home = "/home/birkiy"
    envR = f"{home}/anaconda3/envs/CisGraph/bin/Rscript"

projectRoot = f"{home}/github/CisGraph/Vers2.0"
dataRoot = f"{home}/github/Data/CisGraph/Vers2.0"
figureRoot = f"{home}/github/Figures/CisGraph/Vers2.0"

import os

os.chdir(home)


import csv
import pickle
import subprocess

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

from _collections import deque
import collections
import bisect as bi
import time
