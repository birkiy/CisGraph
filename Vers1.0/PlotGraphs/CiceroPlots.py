

from Functions.PlotFunctions import *
from BFS.CiceroBFS import *



setDhtRoot = set(dhtRooted.keys())
setEthRoot = set(ethRooted.keys())

# Both
both = setDhtRoot.intersection(setEthRoot)

for root in both:
    fig = plt.figure(figsize=[15,12])
    gs = gridspec.GridSpec(ncols=2, nrows=1)

    axDht = fig.add_subplot(gs[0])
    gl = dhtRooted[root]
    sl = dhtShell[root]
    SingleGraphPlot(G, axDht, root, gl, sl)

    axEth = fig.add_subplot(gs[0])
    gl = ethRooted[root]
    sl = ethShell[root]
    SingleGraphPlot(G, axEth, root, gl, sl)

    fig.savefig(f"{figureRoot}/Cicero/{root}.pdf"
