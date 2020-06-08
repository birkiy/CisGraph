
from Functions.CustomBFS import *
from PlotGraphs.PlotRoot import root, depthLimit



# G = ethG
tTree = {}

Shell = {}
Rooted = {}


shell = [None for _ in range(depthLimit+1)]
currentDepth = []

g = bfs_treeCustom(G, root, depth_limit=depthLimit, current_depth=currentDepth, distance=False,shell=shell)


lvlN = depthLimit - min(currentDepth)
if not lvlN in tTree:
    tTree[lvlN] = {}
    if not root in tTree[lvlN]:
        tTree[lvlN][root] = [g]
    else:
        tTree[lvlN][root] += [g]
else:
    if not root in tTree[lvlN]:
        tTree[lvlN][root] = [g]
    else:
        tTree[lvlN][root] += [g]


shell = [i for i in shell if i]

Shell[root] = shell
Rooted[root] = list(g.nodes())
