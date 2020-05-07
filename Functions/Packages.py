

colorPalette = {"arp": "#F99FA5",
                "ctc": "#9FE2D4",
                "fox": "#C1A2EF",
                "med": "#9ECCF2",
                "enz": "#F7C29C",
                "ari": "#F9E0A1",
                "pol": "#F29ECC",
                "nmy": "#E2B79D",

                "gen": "#ddd8c4",

                "upP": "#95342c",
                "dwP": "#2c6589",


                "enP": "#F4A18C",

                "enU": "#fb8182",
                "enD": "#668698",


                "enh": "#a3c9a8",

                "con": "#FADE89",
                "ind": "#57A4B1",
                "non": "#B0D894",


                "tad": "#84b59f",
                "com": "#69a297" ,
                "chr": "#50808e"
}

server = False

if server:
    home = "/kuacc/users/ualtintas20/CisGraph"
    env = "/kuacc/users/ualtintas20/tools/anaconda3/envs/CisGraph"

else:
    home = "/home/birkiy/github/CisGraph"
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
from pygraphviz import *

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
