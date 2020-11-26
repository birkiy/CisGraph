

import pandas as pd

import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import ticker as mticker
from statannot import add_stat_annotation


header = ["chr", "start", "end", "identifier", "mean_signal", "numsamples", "index.ARBS"]
con = pd.read_table("/groups/lackgrp/ll_members/berkay/DHS/creOverARBS/con.cre.bed", names=header)
con["nodeClass"] = "con"
ind = pd.read_table("/groups/lackgrp/ll_members/berkay/DHS/creOverARBS/ind.cre.bed", names=header)
ind["nodeClass"] = "ind"
non = pd.read_table("/groups/lackgrp/ll_members/berkay/DHS/creOverARBS/non.cre.bed", names=header)
non["nodeClass"] = "non"

cre = pd.concat([con, ind, non])
cre["nlog"] = np.log(cre["numsamples"])



fig = plt.figure(figsize=[5,5])
# gs = gridspec.GridSpec(ncols=2, nrows=2)
# plt.subplots_adjust(wspace=0.4, hspace=0.4)

# params = dict(
#     data=tGR,
#     hue="nodeClass",
#     s=8
# )

boxPairs = [("con", "ind"),
            ("con", "non"),
            ("ind", "non")
            ]
ax = sns.violinplot(data=cre, y="nlog", x="nodeClass", palette = ["#63b7af", "#abf0e9", "#d4f3ef"], cut=0)
# plt.yscale('log')
plt.xlabel("NodeClasses")
plt.ylabel("# samples CRE found in")
ax.yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
ax.yaxis.set_ticks([np.log10(x) for p in range(0,10) for x in np.linspace(10**p, 10**(p+1), 10)], minor=True)
add_stat_annotation(ax, data=cre, y="nlog", x="nodeClass", box_pairs=boxPairs, test='Mann-Whitney')


fig.savefig("/groups/lackgrp/ll_members/berkay/DHS/creOverARBS/creARBS.pdf")
