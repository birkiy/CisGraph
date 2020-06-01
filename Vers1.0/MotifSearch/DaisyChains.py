
from Functions.Helpers import *

from Functions.Packages import *

home = "/home/birkiy/github/CisGraph"

C = pickle.load(open(f"{dataRoot}/tmpData/GraphsCData.p", "rb" ))
T = pickle.load(open(f"{dataRoot}/tmpData/GraphsTData.p", "rb" ))
G = pickle.load(open(f"{dataRoot}/tmpData/GraphsGData.p", "rb" ))
P = pickle.load(open(f"{dataRoot}/tmpData/GraphsPData.p", "rb" ))


upP = [_[0] for _ in G.nodes(data="nodeClass") if _[1] in ("enU", "upP")]



source = "FKBP5"

rangeX = G.nodes[source]["nodeRange"]
Gbed = dict(G.nodes(data="nodeRange"))
Gbed = {k: v for k,v in sorted(Gbed.items(), key=lambda kv: kv[1][1])}

gl = rangesFromUpperRange(rangeX[0], rangeX[1] - 200*1000, rangeX[2] + 200*1000,[Gbed],["G"] )
gl = [_[0] for _ in gl if not G.nodes[_[0]]["nodeClass"] in ("gen", "upP", "dwP")]
gl.append(source)


g = G.subgraph(gl).copy()


allPaths = []
for node in g.nodes():

    paths = list(nx.all_simple_paths(g,source, node))
    allPaths += paths




allPaths2 = allPaths.copy()


for i, path1 in enumerate(allPaths):
    for path2 in allPaths[i+1:]:


        print(path2, type(path2))

        tmp1 = set(path1)
        tmp2 = set(path2)

        if tmp2.issubset(tmp1):
            try:
                allPaths2.remove(path2)
            except:
                pass

allPaths3 = allPaths2.copy()

for i, path1 in enumerate(allPaths2[::-1]):
    for path2 in allPaths2[::-1][i+1:]:


        tmp1 = set(path1)
        tmp2 = set(path2)

        if tmp2.issubset(tmp1):
            try:
                allPaths3.remove(path2)
            except:
                pass

allPaths2 = allPaths3
print(len(allPaths2), len(allPaths))

allPaths3 = []
for path in allPaths2:
    starts = [G.nodes[node]["nodeRange"][1] for node in path]
    order = [x for _,x in sorted(zip(starts,path))]

    allPaths3 += [order]

allPaths2 = allPaths3
