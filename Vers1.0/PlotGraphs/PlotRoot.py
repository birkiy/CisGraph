


from Functions.PlotFunctions import *

root = "STEAP4"
depthLimit = 3
from BFS.SingleNodeBFS import *




fig = plt.figure(figsize=[15,15])
gs = gridspec.GridSpec(ncols=1, nrows=1)
#
ax = fig.add_subplot(gs[0])
gl = Rooted[root]
sl = Shell[root]
SingleGraphPlot(G, ax, root, gl, sl)


fig.savefig(f"{figureRoot}/{root}.{depthLimit+1}lvl.BFS.pdf")
plt.close("all")
