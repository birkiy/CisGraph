


from Functions.CustomBFS import *

ethG = pickle.load(open(f"{dataRoot}/tmpData/GraphsG.EtOH.Data.p", "rb" ))
dhtG = pickle.load(open(f"{dataRoot}/tmpData/GraphsG.DHT.Data.p", "rb" ))

Gs = [dhtG, ethG]

# for G in Gs:
G = Gs[0]
tTree = {}
dhtShell = {}
dhtRooted = {}
depthLimit = 4
pro = [_ for _ in G.nodes() if G.nodes[_]["nodeClass"] in ("upP", "dwP")]
for root in pro:
    shell = [None for _ in range(depthLimit+1)]
    currentDepth = []
    g = bfs_treeCustom(G, root, depth_limit=depthLimit, current_depth=currentDepth, distance=False,shell=shell)
    # if len(g) == 1:
    #     continue
    try:
        lvlN = depthLimit - min(currentDepth)
    except:
        lvlN = depthLimit
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
    dhtShell[root] = shell
    dhtRooted[root] = list(g.nodes())


# for G in Gs:
G = Gs[1]
tTree = {}
ethShell = {}
ethRooted = {}
depthLimit = 4
pro = [_ for _ in G.nodes() if G.nodes[_]["nodeClass"] in ("upP", "dwP")]
for root in pro:
    shell = [None for _ in range(depthLimit+1)]
    currentDepth = []
    g = bfs_treeCustom(G, root, depth_limit=depthLimit, current_depth=currentDepth, distance=False,shell=shell)
    # if len(g) == 1:
    #     continue
    try:
        lvlN = depthLimit - min(currentDepth)
    except:
        lvlN = depthLimit
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
    ethShell[root] = shell
    ethRooted[root] = list(g.nodes())
