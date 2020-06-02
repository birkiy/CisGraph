

from Functions.PlotFunctions import *
from BFS.CiceroBFS import *



setDhtRoot = set(dhtRooted.keys())
setEthRoot = set(ethRooted.keys())

# Both
both = setDhtRoot.intersection(setEthRoot)

for root in both:
    fig = plt.figure(figsize=[15,12])
    gs = gridspec.GridSpec(ncols=2, nrows=1)
    #
    axDht = fig.add_subplot(gs[0])
    gl = dhtRooted[root]
    sl = dhtShell[root]
    SingleGraphPlot(Gs[0], axDht, root, gl, sl)
    #
    axEth = fig.add_subplot(gs[1])
    gl = ethRooted[root]
    sl = ethShell[root]
    SingleGraphPlot(Gs[1], axEth, root, gl, sl)
    #
    fig.savefig(f"{figureRoot}/Cicero/{root}.pdf")
    plt.close("all")


# Here, all promoter centered components (left DHT, right EtOH), there are some also duplicates (for example, FKBP5, ARMC12 and LOC285847 are in the same component). 
