

colorPalette = {"arp": "#F99FA5",
                "ctc": "#9FE2D4",
                "fox": "#C1A2EF",
                "med": "#9ECCF2",
                "enz": "#F7C29C",
                "ari": "#F9E0A1",
                "pol": "#F29ECC",
                "nmy": "#E2B79D",

                "pro": "#000000",

                "upP": "#95342c",
                "dwP": "#2c6589",


                # "enP": "#F4A18C",
                #
                # "enU": "#fb8182",
                # "enD": "#668698",


                # "enh": "#a3c9a8",

                "con": "#5A5A5A",
                "ind": "#F9746D",
                "non": "#ACACAC",


                "tad": "#84b59f",
                "com": "#69a297" ,
                "chr": "#50808e",

                # "cre": "#7ebdb4"
}


conditionCheck = {"createGIs": False,

                  "changeGIs": False,

                  "levelInit": True,
                  "levelPlugIN": False,
                  "levelConnect": False,

                  "plotTadBFS": False,

                  "motifSource": False


}


server = False

if server:
    home = "/kuacc/users/ualtintas20/CisGraph/Vers1.0"
    env = "/kuacc/users/ualtintas20/tools/anaconda3/envs/CisGraph"

else:
    home = "/home/birkiy/github/CisGraph/Vers1.0"
    env = "/home/birkiy/anaconda3/envs/CisGraph"

import csv
import pickle
import subprocess


import numpy as np
import networkx as nx
import pandas as pd

from matplotlib import pyplot, patches
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import gridspec, colors
# from pygraphviz import *
from statannot import add_stat_annotation

import matplotlib
import itertools
import scipy
import math
import random
import heapq
from collections import deque




from _collections import deque
import collections
import bisect as bi
